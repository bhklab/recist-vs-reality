import pandas as pd 
import SimpleITK as sitk
import numpy as np
import itertools

from damply import dirs 
from pathlib import Path

def find_multi_ann(matched_data: pd.DataFrame): 
    '''
    Find all instances of a tumour having multiple annotators from the matched data information. 

    Parameters
    ----------
    matched_data: pd.DataFrame 
        Contains the matched RECIST measurements (annotations), images and segmentations and their associated information 

    Returns
    ----------
    mult_ann_data: pd.DataFrame
        A subset of the matched data that only has tumours with multiple annotations.
    ''' 
    # Find all duplicates within the SegNIFTILocation column (if two annotations map to the same segmentation, this is an instance of multiple annotators)
    mult_ann_data = matched_data[matched_data.duplicated(subset = ['SegNIFTILocation'], keep = False)]

    return mult_ann_data

def get_vox_vol(rel_nifti_path: Path, 
                dataset: str, 
                disease_site: str): 
    '''
    Gets the voxel volume in cubic mm from a given NIFTI segmentation file. 

    Parameters
    ----------
    rel_nifti_path: Path 
        Relative path to the NIFTI segmentation file starting from the disease site.

    Returns 
    ----------
    vox_vol: int 
        The voxel volume (in cubic mm)
    '''
    full_nifti_path = dirs.PROCDATA / disease_site / dataset / 'images' / rel_nifti_path

    # Load in NIFTI file, get voxel spacing, and convert to array
    nifti_img = sitk.ReadImage(full_nifti_path)
    nifti_arr = sitk.GetArrayFromImage(nifti_img)
    vox_spacing = nifti_img.GetSpacing()
    vox_unitary_vol = vox_spacing[0] * vox_spacing[1] * vox_spacing[2]

    # Count all of the voxels in the segmentation 
    vox_count = np.count_nonzero(nifti_arr) 
    
    # Calculate voxel volume 
    vox_vol = vox_count * vox_unitary_vol 

    return vox_vol

def comp_pairwise_info(multiple_ann_df: pd.DataFrame,
                       dataset: str, 
                       disease_site: str): 
    '''
    Compute all pairwise combinations (at the individual tumour level) for the multiple annotations and save them into a dataframe.

    Parameters
    ----------
    multiple_ann_df: pd.DataFrame 
        A dataframe only containing the tumours which had multiple annotations. 
    dataset: str
    disease_site: str

    Returns
    ----------
    mult_ann_calcs_df: pd.DataFrame 
        Contains calculations made pairwise between different annotations and important information about the segmentations. 
    '''
    cols = ['Ann1Length', 
            'Ann2Length', 
            'AnnLengthDiff', 
            'AnnAbsLengthDiff', 
            'AnnLengthRatio', 
            'AnnLengthPercDiff', 
            'AnnAbsLengthPercDiff',
            'SegVoxVolume',
            'meta_PatientID',
            'meta_SegNIFTIShortLoc', 
            'meta_Ann1SeriesUID',
            'meta_Ann2SeriesUID',
            'meta_Ann1RefSlice', 
            'meta_Ann2RefSlice']
    
    # Get a list of all duplicate values for SegNIFTILocation to subset and iterate over
    unique_seg_vals = multiple_ann_df['SegNIFTILocation'].values.tolist()

    for seg in unique_seg_vals: 
        curr_tum_set = multiple_ann_df[multiple_ann_df['SegNIFTILocation'] == seg] #Subset to lesion specific multiple annotators

        # Get volume for the current segmentation 
        curr_seg_rel_path = curr_tum_set['SegNIFTIShortLoc'].iloc[0] #These should all be the same so I can take the first instance of it 
        curr_voxel_volume = get_vox_vol(rel_nifti_path = curr_seg_rel_path, 
                                        dataset = dataset, 
                                        disease_site = disease_site)

        # Get all unique combinations of multiple annotations 
        num_rows = curr_tum_set.shape[0]
        iloc_list = list(range(num_rows))
        unique_pairs = list(itertools.combinations(iloc_list, 2))

        # Iterate over all unique combinations of annotations 
        for pair in unique_pairs: 
            ann_1_info = curr_tum_set.iloc[pair[0]]
            ann_2_info = curr_tum_set.iloc[pair[1]]
            
            ann_1_length = ann_1_info['AnnLongAxisLength']
            ann_2_length = ann_2_info['AnnLongAxisLength']

            # Calculate the RECIST measurement 
            ann_length_diff = ann_2_length - ann_1_length
            ann_abs_length_diff = abs(ann_length_diff)
            ann_length_ratio = min(ann_1_length, ann_2_length) / max(ann_1_length, ann_2_length)

            ann_length_perc_diff = (ann_2_length - ann_1_length) / ((ann_1_length + ann_2_length)/2) * 100
            ann_abs_len_perc_diff = abs(ann_length_perc_diff)

            # Get all relevent metadata from both 
            patient_id = ann_1_info['PatientID']
            ann1_serUID = ann_1_info['AnnSeriesInstanceUID']
            ann2_serUID = ann_2_info['AnnSeriesInstanceUID']
            ann1_sliceID = ann_1_info['AnnReferencedSOPUID']
            ann2_sliceID = ann_2_info['AnnReferencedSOPUID']
            curr_info = [ann_1_length, 
                         ann_2_length, 
                         ann_length_diff, 
                         ann_abs_length_diff,
                         ann_length_ratio, 
                         ann_length_perc_diff, 
                         ann_abs_len_perc_diff, 
                         curr_voxel_volume,
                         patient_id, 
                         seg,
                         ann1_serUID, 
                         ann2_serUID,
                         ann1_sliceID,
                         ann2_sliceID]
            
            if 'mult_ann_calcs_df' not in locals(): 
                mult_ann_calcs_df = pd.DataFrame([curr_info], columns = cols, index = [0])
            else: 
                temp = pd.DataFrame([curr_info], columns = cols, index = [0])
                mult_ann_calcs_df = pd.concat([mult_ann_calcs_df, temp], ignore_index = True).reset_index(drop = True)

    return mult_ann_calcs_df

def run_multi_ann(dataset, 
                  disease_site): 
    matched_path = dirs.PROCDATA / disease_site / dataset / Path('metadata/annotation_seg_matching/matching_ann_to_seg.csv')
    matched_data = pd.read_csv(matched_path)
    multi_df = find_multi_ann(matched_data = matched_data)

    multi_calc_df = comp_pairwise_info(multiple_ann_df = multi_df, 
                                       dataset = dataset, 
                                       disease_site = disease_site) 
    #Drop dupes 
    multi_calc_df = multi_calc_df.drop_duplicates()
    
    multi_calc_df.to_csv('test_multi_ann_calc_NSCLC.csv')

if __name__ == "__main__": 
    run_multi_ann(dataset = 'TCIA_NSCLC-Radiogenomics', 
                  disease_site = 'Lung')
    