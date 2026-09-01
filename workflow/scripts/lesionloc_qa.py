# Quick script to filter results of lesion locator segmentations based on which patients started with a mask that had less than 10 voxels

import pandas as pd 
import numpy as np 
import click
import SimpleITK as sitk

from pathlib import Path

def find_drop_segs(aaura_idx_path: Path, 
                   vox_vol_thresh: int): 
    '''
    Filters out segmentations that have less than the voxel volume threshold provided and returns them in a list.

    Parameters
    ----------
    aaura_idx_path: Path
        The path to the AAuRA index file. Must contain the headers "mask_path" and "mask_volume". Must be csv
    vox_vol_thresh: int
        The voxel volume threshold for filtering. Default is 10. 

    Returns
    ----------
    seg_to_drop: list 
        A list of mask paths that have a voxel volume below the threshold given. 
    '''
    # Load in csv file 
    aaura_idx = pd.read_csv(aaura_idx_path)

    # Subset samples that do not meet the voxel volume threshold
    seg_to_drop_df = aaura_idx[aaura_idx['mask_volume'] < vox_vol_thresh]

    # Save out dropped segmentations in the same location as the AAuRA index 
    dropped_savepath = Path("/".join(aaura_idx_path.split("/")[:-1])) / Path('dropped_segmentations_' + str(vox_vol_thresh) + 'voxvol.csv')
    seg_to_drop_df.to_csv(dropped_savepath, index = False)

    # Get mask paths that achieve that threshold and return 
    seg_to_drop = seg_to_drop_df['mask_path'].values.tolist()

    return seg_to_drop 

def list_nonzero_seg_slices(seg: np.ndarray): 
    '''  
    From a given 3D segmentation array, list the slices that have nonzero values (mask) in them.

    Parameters
    ----------
    seg: np.ndarray
        A 3D array containing a mask (ground truth, predicted, etc.)
    
    Returns
    ----------
    nonzero_slices: list 
        Contains all of the slice numbers where there are nonzero values
    '''
    nonzero_slices = []
    for slice_idx in range(seg.shape[0]): 
        if np.count_nonzero(seg[slice_idx]) > 0: 
            nonzero_slices.append(slice_idx)
    return nonzero_slices

def subset_results_slices(mask_subset: list, 
                          results_path: Path): 
    ''' 
    From a list of mask paths, check to see which patients need to be dropped based on matching the slices found in the segmentation. 

    Parameters 
    ----------
    mask_subset: list 
        The list of mask paths to drop 
    results_path: Path
        The path to the inference metrics results. Must be csv
    
    Returns
    ---------
    results_df: pd.DataFrame
        The new dataframe with all specified segmentations dropped
    '''
    # Load in the results csv 
    results_df = pd.read_csv(results_path) 

    # Go through each of the paths in the list, find mask slice range, and compare to the results dataframe based on filename and slice range
    for mask_path in mask_subset: 
        curr_mask = sitk.ReadImage("data/procdata/MultiSite/" + mask_path)
        curr_mask_arr = sitk.GetArrayFromImage(curr_mask)

        gts_slices = str(list_nonzero_seg_slices(curr_mask_arr))

        # Get current Patient ID and row to drop
        patient_ID = mask_path.split("/")[-2]
        row_to_drop = results_df[(results_df['filename'].str.contains(patient_ID, case=True)) & (results_df['GTSliceList'] == gts_slices)]
        
        # Drop the row
        results_df = results_df.drop([row_to_drop.index.tolist()[0]])

    return results_df

def align_gt_pred(filtered_idx: pd.DataFrame, 
                  filtered_results: pd.DataFrame): 
    """  
    From the filtered AAuRA index and prediction results file, align the ground truth and predicted 
    results into one dataframe.

    Parameters
    ----------
    filtered_idx: pd.DataFrame
        The AAuRA index file after it has been filtered by the lower bound threshold 
    filtered_results: pd.DataFrame
        The results after the corresponding small volume masks have been dropped

    Returns
    ----------
    gt_pred_df: pd.DataFrame
        The results and ground truth data aligned to be one row for each mask
    """

    #Iterate over all masks within the index file
    for mask_path in filtered_idx['mask_path'].values.tolist(): 
        # Get current row from the index and the current mask information from the NIfTI file 
        curr_idx_row = filtered_idx[filtered_idx['mask_path'] == mask_path]

        curr_mask = sitk.ReadImage('data/procdata/MultiSite/' + mask_path) 
        curr_mask_arr = sitk.GetArrayFromImage(curr_mask) 

        # Get current Patient ID and row to align
        path_match = "/".join(mask_path.split("/")[-3:]).removesuffix('.nii.gz') + "_pred.nii.gz" 
        row_to_align = filtered_results[filtered_results['filename'].str.contains(path_match, case=True)]

        # Check if there is no match, make note of it (this shouldn't happen but just in case)
        if row_to_align.empty: 
            print(f"Mask path {mask_path} does not have a match in the predicted data. Please debug manually.")
            continue

        # Drop that row from the copy of the results to keep track of samples that may have been dropped from ground truth but still ran 
        filtered_results = filtered_results.drop([row_to_align.index.tolist()[0]]) 

        # Join the two rows horizontally 
        temp_df = pd.concat([curr_idx_row.reset_index(drop = True), row_to_align.reset_index(drop = True)], axis = 1)
        if 'gt_pred_df' not in locals():
            gt_pred_df = temp_df
        else:
            gt_pred_df = pd.concat([gt_pred_df, temp_df], ignore_index = True)
            print(f"Appended {mask_path} info")
    if not filtered_results.empty: 
        print("Some predicted segmentations now have no ground truth because of the previous filtering.") 
        filtered_results.to_csv("preds_w_missing_gt.csv") 

    gt_pred_df = gt_pred_df.reset_index(drop = True) 

    return gt_pred_df

@click.command() 
@click.option('--aaura_idx_path')
@click.option('--vox_vol_thresh', default = 10)
@click.option('--results_path') 
@click.option('--align_idx_results', default = True)
def filter_results(aaura_idx_path: Path, 
                   vox_vol_thresh: int, 
                   results_path: Path, 
                   align_idx_results: bool): 
    '''  
    From the AAuRA index file, find all segmentations with a voxel volume below the specified threshold and remove it from the results file. 

    Parameters
    ----------
    aaura_idx_path: Path
        The path to the AAuRA index file. Must contain the headers "mask_path" and "mask_volume". Must be csv
    vox_vol_thresh: int
        The voxel volume threshold for filtering.
    results_path: Path
        The path to the inference metrics results. Must be csv
    align_idx_results: bool
        Combine the ground truth (idx_path) and the predicted (results) into one dataframe 
    '''
    drop_seg_list = find_drop_segs(aaura_idx_path = aaura_idx_path, 
                                   vox_vol_thresh = vox_vol_thresh)
    
    filt_results = subset_results_slices(mask_subset = drop_seg_list, 
                                         results_path = results_path) 

    # If align is True, align and export the data as one csv 
    if align_idx_results:
        idx_dropped_savepath = Path("/".join(aaura_idx_path.split("/")[:-1])) / Path('dropped_segmentations_' + str(vox_vol_thresh) + 'voxvol.csv')
        drop_idx = pd.read_csv(idx_dropped_savepath)

        aaura_idx_df = pd.read_csv(aaura_idx_path)
        idx_keep_rows = ~aaura_idx_df['mask_path'].isin(drop_idx['mask_path'])

        filtered_idx = aaura_idx_df[idx_keep_rows].reset_index(drop = True) 

        aligned_gt_pred_df = align_gt_pred(filtered_idx = filtered_idx, 
                                           filtered_results = filt_results)
        
        aligned_gt_pred_path = Path(str("/".join(results_path.split("/")[:-1]))) / Path('gt_pred_drop' + str(vox_vol_thresh) + 'voxvol.csv')
        aligned_gt_pred_df.to_csv(aligned_gt_pred_path, index = False)

    # Reset the index of the results dataframe and save 
    filt_results = filt_results.reset_index(drop=True) 
    filt_result_path = Path(str("/".join(results_path.split("/")[:-1]))) / Path('metric_eval_drop' + str(vox_vol_thresh) + 'voxvol.csv')
    filt_results.to_csv(filt_result_path, index = False)

if __name__ == "__main__": 
    filter_results()
