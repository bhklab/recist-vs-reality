""" Script to get RECIST measurement changes from the baseline to synthetic masks
and generate the visualizations."""

import pandas as pd
import numpy as np
import SimpleITK as sitk
import click

from imgtools.transforms.functional import resample
from tqdm import tqdm
from joblib import Parallel, delayed
from pathlib import Path
from skimage.measure import regionprops, label
from skimage.draw import line

def get_vox_vol(mask_array: np.ndarray, unitary_vol: float):
    """
    Gets the voxel volume in cubic mm from a given NIFTI segmentation file.

    Parameters
    ----------
    mask_array: np.ndarray
        The array of the 3D binary mask of a tumour
    unitary_vol: float
        The volume of a single voxel for this specific mask array

    Returns
    ----------
    vox_vol: int
        The voxel volume (in cubic mm)
    """
    # Count all of the voxels in the segmentation
    vox_count = np.count_nonzero(mask_array)

    # Calculate voxel volume
    vox_vol = vox_count * unitary_vol

    return vox_vol

def pad_bbox(box:np.array,
             mask:np.array,
             padding:int,
             spacing:np.array = None
             ) -> np.array:
    # Get full image dimensions to keep padding within image size
    # D, H, W (z, y, x)
    mask_shape = mask.shape

    if spacing is not None: # Use the actual image spacing to calculate the padding
        # Check that spacing can be applied to this mask's bounding box dimensions
        if len(spacing) < len(mask_shape):
            spacing = spacing[0, len(mask_shape)]
        elif len(spacing) > len(mask_shape):
            message = "Spacing for padding has more dimensions than the mask image."
            raise ValueError(message)

        # calculate the number of voxels to pad based on the actual image spacing
        padding = np.round(padding / spacing)
    else:
        # Convert padding into an array with the same length as image dimensions
        # This matches the behaviour of the spacing option
        padding = padding * np.ones(len(mask_shape))

    pad_x_min = max(0, box[0] - padding[0])
    pad_y_min = max(0, box[1] - padding[1])
    # Handling 2D bounding box
    if len(box) == 4:
        mask_H, mask_W = mask_shape[0], mask_shape[1]
        pad_x_max = min(mask_W, box[2] + padding[0])
        pad_y_max = min(mask_H, box[3] + padding[1])

        padded_box = np.array([pad_x_min, pad_y_min,
                               pad_x_max, pad_y_max])

    # Handling 3D bounding box
    if len(box) == 6:
        mask_D, mask_H, mask_W = mask_shape[0], mask_shape[1], mask_shape[2]
        pad_z_min = max(0, box[2] - padding[2])
        pad_x_max = min(mask_W, box[3] + padding[0])
        pad_y_max = min(mask_H, box[4] + padding[1])
        pad_z_max = min(mask_D, box[5] + padding[2])

        padded_box = np.array([pad_x_min, pad_y_min, pad_z_min,
                               pad_x_max, pad_y_max, pad_z_max])

    return padded_box.astype(int)

def mask2D_to_bbox(gt2D:np.array,
                   mask_path:Path,
                   padding:int | None = None,
                   spacing:np.array = None
                   ) -> np.array:
    try:
        ## Old code
        # y_indices, x_indices = np.where(gt2D > 0)
        # x_min, x_max = np.min(x_indices), np.max(x_indices)
        # y_min, y_max = np.min(y_indices), np.max(y_indices)
        # boxes = np.array([x_min, y_min, x_max, y_max])

        props = regionprops(gt2D)[0]
        y_cent, x_cent = props.centroid
        orientation = props.orientation
        semi_maj_axis_len = props.major_axis_length / 2

        x_start = x_cent - np.sin(orientation) * semi_maj_axis_len
        y_start = y_cent - np.cos(orientation) * semi_maj_axis_len

        x_end = x_cent + np.sin(orientation) * semi_maj_axis_len
        y_end = y_cent + np.cos(orientation) * semi_maj_axis_len

        boxes = np.array([x_start, y_start, x_end, y_end])

        if padding:
            boxes = pad_bbox(box = boxes,
                             mask = gt2D,
                             padding = padding,
                             spacing = spacing)

        return boxes.astype(int)

    except Exception as e:
        raise Exception(f'error {e} with file {mask_path} and sum of gts is {gt2D.sum()}')

def get_line_from_recist(recist_coords: np.array,
                         slice_number: int,
                         img_size: np.array):
    '''
    From the RECIST measurement coordinates, generate a line connecting both coordinates on the correct slice and return an np.ndarray the same shape as the image.
    Output to be compatible with the ['recist'] array of the .npz files needed for MedSAM2-RECIST.

    Parameters
    ----------
    recist_coords: array
        A list of coordinates in [x1, y1, x2, y2] format that defines the RECIST measurement
    slice_number: int
        The slice that the measurement was taken on
    img_size: np.array
        The x, y, and z size of the image in [z_space, x_space, y_space] format

    Returns
    ----------
    recist_arr: np.ndarray
        A binary array of the same shape as the image with the pixels of the line = 1
    '''
    # Check to see if RECIST coordinates are in string form and if so, convert to list
    if type(recist_coords) == str:
        just_coords = recist_coords.strip("[]")
        recist_coords = list(map(float, just_coords.split()))
        print(recist_coords)

    #Generate an array in the same size as the image filled with all zeros
    recist_arr = np.zeros((img_size[0], img_size[1], img_size[2]), dtype = int)

    #Round the coordinate values to their nearest integers
    coords_round = np.rint(recist_coords).astype(int)

    #Draw line using coordinates
    rr, cc = line(coords_round[0], coords_round[1], coords_round[2], coords_round[3])

    #Put line into the correct slice in the RECIST array of all zeros
    recist_arr[slice_number][cc, rr] = 1

    return recist_arr

def get_recist_tum(diam_change: float):
    '''
    Get the RECIST category from the diameter change.

    Parameters
    ----------
    diam_change: float
        The percent diameter change of from the first to second timepoint masks

    Returns
    ----------
    recist_cat: str
        Either PD, SD, PR, or CR
    '''
    if diam_change == None: 
        return None
    elif diam_change == -100: 
        return 'CR'
    elif -100 < diam_change <= -30:
        return 'PR'
    elif diam_change > -30 and diam_change < 20:
        return 'SD'
    elif diam_change > 20:
        return 'PD'
    else: 
        return None

def get_vol_recist(vol_change: float): 
    '''  
    Get the volumetric RECIST category (as defined here: https://pubs.rsna.org/doi/full/10.1148/rycan.220166#:~:text=Categorical%20response%20thresholds%20for%20volumetric,of%20volumes%20(Fig%20S1)) from the volume change. 

    Parameters
    ----------
    vol_change: float 
        The percentage volume change from the first to second timepoint masks 
    
    Returns 
    ----------
    vol_recist_at: str
        Either PD, SD, PR, or CR
    '''
    if vol_change == None: 
        return None
    elif vol_change == -100: 
        return 'CR'
    elif -100 < vol_change <= -65: 
        return 'PR'
    elif vol_change > -65 and vol_change <= 73: 
        return 'SD' 
    elif vol_change > 73: 
        return 'PD'
    else: 
        return None

def calc_recist_info(mask_path: Path): 
    """  
    From a mask image, calculate the RECIST diameter and the voxel volume based on the slice 
    with the largest pixel area. 

    Parameters
    ----------
    mask_path: Path
        The path to the current mask 
    
    Returns
    ---------- 
    recist_len: float
        The largest axial diameter found 
    voxel_vol: float 
        The voxel volume of the current mask 
    """
    # Read in mask image
    mask_img = sitk.ReadImage(mask_path) 

    # Get the spacing from the image for further calculations 
    spacing = mask_img.GetSpacing() 

    # Find the slice with the largest area for further measurement 
    mask_arr = sitk.GetArrayFromImage(mask_img) 
    
    slice_sums = np.sum(mask_arr, axis = (1,2))
    max_area_slice = np.argmax(slice_sums) 

    # Check if the mask is empty, if it is, return zeros for the length and volume 
    if mask_arr[max_area_slice].sum() == 0: 
        return 0, 0
    
    # Get the pixel RECIST coordinates of the max slice 
    bbox_2d = mask2D_to_bbox(gt2D = mask_arr[max_area_slice], 
                               mask_path = 'not_applicable')
    
    # Convert pixel coordinates to physical coordinates 
    bbox_coord1 = (float(bbox_2d[0]), float(bbox_2d[1]), float(max_area_slice)) 
    bbox_coord2 = (float(bbox_2d[2]), float(bbox_2d[3]), float(max_area_slice))

    print(bbox_coord1) 
    print(bbox_coord2) 

    bbox_phys1 = mask_img.TransformContinuousIndexToPhysicalPoint(bbox_coord1) 
    bbox_phys2 = mask_img.TransformContinuousIndexToPhysicalPoint(bbox_coord2)

    # Calculate RECIST measurement length using standard Euclidean distance 
    recist_len = np.linalg.norm(np.array(bbox_phys2) - np.array(bbox_phys1))

    # Calculate the voxel volume of the mask 
    unit_vol = spacing[0] * spacing[1] * spacing[2]
    voxel_vol = get_vox_vol(mask_array = mask_arr, 
                            unitary_vol = unit_vol) 
    
    return recist_len, voxel_vol

def get_indiv_percdiam_change(measure_t1: float, 
                              measure_t2: float): 
    """  
    From the measurements given, calculate the percentage change for a singular tumour. 

    Parameters
    ----------
    measure_t1: float 
        The RECIST/volume measurement length of the first timepoint (in mm or mm^3) 
    measure_t2: float 
        The RECIST/volume measurement length of the second timepoint (in mm or mm^3) 

    Returns
    ----------
    perc_change: float 
        The percentage change of the diameters across the two timepoints
    """
    # Calculate percentage change
    if measure_t1 == 0 and measure_t2 == 0: 
        return None # No segmentation detected here so this wouldn't even get a RECIST criteria and should be dropped later
    elif measure_t1 == 0 and measure_t2 != 0: 
        return 9999999.99 # This simulates a new tumour showing up, should be PD
    else:
        perc_diam_change = (measure_t2 - measure_t1) / measure_t1 * 100 
        return perc_diam_change 

def calc_2d_dice(mask_gt: Path, 
                 mask_pred: Path):
    """  
    From the predicted and ground truth masks, calculate the 2D Dice Score using the largest pixel area slice from the ground truth. 

    Parameters 
    ----------
    mask_gt: Path 
        The path to the ground truth mask. 
    mask_pred: Path 
        The path to the predicted mask. 

    Returns
    ----------
    dice_2d: float
        The 2D Dice measurement calculated. 
    """
    # Read in both masks and convert to arrays
    gt_img = sitk.ReadImage(mask_gt) 
    pred_img = sitk.ReadImage(mask_pred) 

    gt_arr = sitk.GetArrayFromImage(gt_img) 
    pred_arr = sitk.GetArrayFromImage(pred_img) 

    # Find the largest slice in the ground truth mask 
    slice_sums = np.sum(gt_arr, axis = (1,2))
    max_area_slice = np.argmax(slice_sums) 

    gt_large_slice = gt_arr[max_area_slice] 
    pred_large_slice = pred_arr[max_area_slice]

    # Calculate and return the 2D Dice Score
    dice_2d = (np.sum(pred_large_slice[gt_large_slice == 1]) * 2.0) / (np.sum(pred_large_slice) + np.sum(gt_large_slice))

    return dice_2d

def recist_eval_indiv(t1_mask_path: Path,
                      t2_mask_path: Path) -> pd.DataFrame: 
    """  
    From a given timepoint one and two mask pair, calculate the voxel volume and largest diameter (in mm) for both and 
    get lesion-level RECIST category along with the percentage volume change. 

    Parameters 
    ----------
    t1_mask_path: Path
        The path to the first timepoint mask 
    t2_mask_path: Path 
        The path to the second timepoint mask 

    Returns
    ----------
    recist_summary: pd.DataFrame
        A dataframe containing information on the mask paths and RECIST/volume change 
        statistics.
    """
    # Calculate the RECIST and volume information for the two timepont masks
    t1_recist_len, t1_vox_vol = calc_recist_info(mask_path = t1_mask_path) 
    t2_recist_len, t2_vox_vol = calc_recist_info(mask_path = t2_mask_path) 

    # Get the percentage change over time for both the RECIST diameter and volume measurements 
    diam_perc_change = get_indiv_percdiam_change(measure_t1 = t1_recist_len, 
                                                 measure_t2 = t2_recist_len) 
    vol_perc_change = get_indiv_percdiam_change(measure_t1 = t1_vox_vol, 
                                                measure_t2 = t2_vox_vol) 
    
    # Calculate the RECIST category for this tumour 
    recist_cat = get_recist_tum(diam_change = diam_perc_change) 

    vol_recist_cat = get_vol_recist(vol_change = vol_perc_change)

    # Summarize all information into a dataframe 
    col_names = ['t1_mask_path',
                 't2_mask_path',
                 't1_max_ax_diam', 
                 't2_max_ax_diam', 
                 't1_vox_vol', 
                 't2_vox_vol', 
                 'perc_diam_change', 
                 'perc_vol_change', 
                 'RECIST_cat',
                 'vol_RECIST_cat'
    ]

    curr_info = [t1_mask_path, 
                 t2_mask_path, 
                 t1_recist_len, 
                 t2_recist_len, 
                 t1_vox_vol, 
                 t2_vox_vol, 
                 diam_perc_change, 
                 vol_perc_change, 
                 recist_cat, 
                 vol_recist_cat]
    
    recist_summary = pd.DataFrame([curr_info], columns = col_names, index = [0]) 
    
    return recist_summary 

def run_one_samp_plot_data(t1_gt_mask_path: Path, 
                           t2_gt_mask_path: Path,
                           t1_pred_mask_path: Path, 
                           t2_pred_mask_path: Path, 
                           t1_DICE: float, 
                           t2_DICE: float, 
                           t1_H95: float, 
                           t2_H95: float, 
                           t1_surf_dist: float, 
                           t2_surf_dist: float) -> pd.DataFrame: 
    """ 
    Get all data necessary for simulation study plotting. 

    Parameters
    ----------
    t1_gt_mask_path: Path
        The path to the ground truth first timepoint mask 
    t2_gt_mask_path: Path
        Path to the ground truth second timepoint mask 
    t1_pred_mask_path: Path 
        The path to the predicted first timepoint segmentation 
    t2_pred_mask_path: Path 
        The path to the predicted second timepoint segmentation 
    t1_DICE: float
        The volumetric DICE performance of the first timepoint segmentation 
    t2_DICE: float 
        The volumetric DICE performance of the second timepoint segmentation 
    t1_H95: float
        The Hausdorff distance for the first timepoint segmentation 
    t2_H95: float 
        The Hausdorff distance for the second timepoint segmentation 
    t1_surf_dist: float
        The surface distance performance of the first timepoint segmentation 
    t2_surf_dist: float
        The surface distance performance of the second timepoint segmentation

    Returns
    ----------
    plotting_data: pd.DataFrame 
        A summary of all information necessary for simulation study plotting, including mask paths, RECIST/volume information, and 
        metric performance
    """
    # Get ground truth statistics 
    gt_recist_summ = recist_eval_indiv(t1_mask_path = t1_gt_mask_path, 
                                       t2_mask_path = t2_gt_mask_path)
    
    gt_recist_summ = gt_recist_summ.add_prefix('GT_')

    # Get predicted statistics 
    pred_recist_summ = recist_eval_indiv(t1_mask_path = t1_pred_mask_path, 
                                         t2_mask_path = t2_pred_mask_path)
    
    pred_recist_summ = pred_recist_summ.add_prefix('PRED_')

    # Calculate 2D Dice Score based for both timepoints 
    t1_dice2d = calc_2d_dice(mask_gt = t1_gt_mask_path, 
                             mask_pred = t1_pred_mask_path)
    
    t2_dice2d = calc_2d_dice(mask_gt = t2_gt_mask_path, 
                             mask_pred = t2_pred_mask_path)

    # Organize the performance statistics into a dataframe 
    cols = ['t1_DICE', 
            't2_DICE', 
            't1_2D_DICE',
            't2_2D_DICE',
            't1_HD95', 
            't2_HD95', 
            't1_surface_distance', 
            't2_surface_distance']
    
    metrics = [t1_DICE, 
               t2_DICE,
               t1_dice2d, 
               t2_dice2d, 
               t1_H95, 
               t2_H95, 
               t1_surf_dist,
               t2_surf_dist]

    metrics_df = pd.DataFrame([metrics], columns = cols, index = [0]) 

    # Concatenate all necessary data together and return 
    plot_data_list = [gt_recist_summ, pred_recist_summ, metrics_df]

    plotting_data = pd.concat(plot_data_list, axis = 1)

    plotting_data['RECIST_flip'] = plotting_data['GT_RECIST_cat'] != plotting_data['PRED_RECIST_cat']

    return plotting_data

def get_merge_data(mask_path): 
    """ 
    Used in order to get a common column to merge on based on the mask paths of t1 and t2 masks 

    Parameters
    ----------
    mask_path: str
        The path to the current mask in question 
    
    Returns
    ----------
    merge_data: str 
        The string that will be used to merge the two aligned dataframes on 
    """
    merge_data = "/".join(mask_path.split("/")[-2:])

    return merge_data 

@click.command()
@click.option('--align_t1_path') 
@click.option('--align_t2_path') 
@click.option('--n_jobs')
def run_get_sim_data(align_t1_path: Path, 
                     align_t2_path: Path, 
                     n_jobs: int): 
    """  
    From the aligned data paths, get all necessary information for simulation study plotting and 
    save into one collective dataframe in results path. 

    Parameters 
    ----------
    align_t1_path: Path 
        The path to the aligned AAuRA index and results file for the first timepoint. Assumes csv. 
    align_t2_path: Path
        The path to the aligned AAuRA index and results file for the second timepoint. Assumes csv.
    n_jobs: int 
        The number of jobs to run in parallel. 
    """
    # Read in the two aligned dataframes 
    align_t1_df = pd.read_csv(align_t1_path) 
    align_t2_df = pd.read_csv(align_t2_path) 

    # Put correct suffix in the filename position
    align_t1_df['filename'] = align_t1_df['filename'].str.replace('.nii.gz', '_pred_PTS_25_75_BBOX_MINAX.nii.gz')
    align_t1_df['filename'] = align_t2_df['filename'].str.replace('.nii.gz', '_pred_PTS_25_75_BBOX_MINAX.nii.gz')

    # Give each of the dataframes a prefix for their columns so they can be concatenated (for parallel processing) 
    # and create a common column to merge on so all information can be aligned correctly
    align_t1_df = align_t1_df.add_prefix('t1_') 
    align_t2_df = align_t2_df.add_prefix('t2_') 

    # If there was no overlap in the segmentations (0 dice), the HD will be blank. Need to fill those blanks with 'inf' to avoid accidental dropping
    align_t1_df['t1_hausdorff'] = align_t1_df['t1_hausdorff'].fillna('inf')
    align_t2_df['t2_hausdorff'] = align_t2_df['t2_hausdorff'].fillna('inf') 

    # Create the merge column 
    align_t1_df['merge_col'] = align_t1_df['t1_mask_path'].apply(get_merge_data)
    align_t2_df['merge_col'] = align_t2_df['t2_mask_path'].apply(get_merge_data)

    aligned_summary = pd.merge(align_t1_df.sort_values(by = 'merge_col'), align_t2_df.sort_values(by = 'merge_col'), on = 'merge_col', how = 'outer')

    aligned_summary.to_csv('output_qa.csv') 
    # Drop any rows that that were not merged properly (i.e. rows with NaN values that did not have matching info in both timepoints)
    aligned_summary = aligned_summary.dropna() 

    #Go through each row in the aligned summary and calculate the necessary information 
    sim_plot_data = Parallel(n_jobs = n_jobs)(delayed(run_one_samp_plot_data)
                                           (t1_gt_mask_path = 'data/procdata/MultiSite/' + row['t1_mask_path'],
                                            t2_gt_mask_path = 'data/procdata/MultiSite/' + row['t2_mask_path'],
                                            t1_pred_mask_path = row['t1_filename'],	
                                            t2_pred_mask_path = row['t2_filename'],
                                            t1_DICE = row['t1_volume_dice'],
                                            t2_DICE = row['t2_volume_dice'],
                                            t1_H95 = row['t1_hausdorff'],
                                            t2_H95 = row['t2_hausdorff'],
                                            t1_surf_dist = row['t1_surface_distance'], 
                                            t2_surf_dist = row['t2_surface_distance'])
                                            for _, row in tqdm(aligned_summary.iterrows(), total = aligned_summary.shape[0]))
    
    for plot_data in sim_plot_data: 
        if 'sim_plot_data_df' not in locals(): 
            sim_plot_data_df = plot_data 
        else: 
            sim_plot_data_df = pd.concat([sim_plot_data_df, plot_data], ignore_index = True)
    
    # Save plotting information 
    sim_plot_data_df.to_csv('ll_sim_plotting_info.csv', index = False)

if __name__ == '__main__': 
    run_get_sim_data()