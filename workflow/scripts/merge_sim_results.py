import pandas as pd 
import numpy as np 
import click
import SimpleITK as sitk 

from damply import dirs 
from pathlib import Path 
from functools import reduce 
from skimage.measure import regionprops, label
from imgtools.transforms.functional import resample 

def common_name_finder(filename: str):
    '''
    For the various filenames in the csvs, will get the common string between them (usually sampleID_rtstructID form) for a specified filename. 
    To be applied to a column.

    Parameters
    ----------
    filename: str
        The name recorded in the csv 

    Returns 
    ----------
    common_name: str
        The common name in the form of sampleID_rtstructID
    '''
    if "/" in filename: 
        filename = filename.split("/")[-1] #The last part of the relative file path will have the actual filename we want to use
    
    # Split filename by "." to get the first element containing the common string pattern 
    common_name = filename.split(".")[0]

    return common_name

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

def locate_centre_slice(mask_3d):
    """
    Locates the center slice of a 3D mask.

    Args:
    - mask_3d (ndarray): A 3D binary mask.

    Returns:
    - (int): The index of the center slice.
    """
    
    # find the slice in the center
    lesion_labels = label(mask_3d)
    centre_slc = regionprops(lesion_labels)[0].centroid[0]

    # number of voxels per slice
    vox_per_slc = np.array([np.sum(mask_3d[slc,:,:]) for slc in range(mask_3d.shape[0])])
    max_vox_slc = np.where(vox_per_slc==np.max(vox_per_slc))[0]
    if len(max_vox_slc) > 1: 
        max_vox_slc = max_vox_slc[int(np.floor(len(max_vox_slc)/2))]

    return int(np.floor((centre_slc+max_vox_slc)/2))

def get_major_axis(mid_slice: np.ndarray): 
	'''
	From the middle slice of a mask, get the major axis length. 

	Parameters
	----------
	mid_slice: np.ndarray
		Binary mask slice 
	
	Returns
	-----------
	maj_axis_len: float 
		The major axis length of the binary mask slice
	'''
	# Get region properties 
	region_info = regionprops(mid_slice.astype(int)) 
	try: 
		maj_axis_len = region_info[0].axis_major_length
	except IndexError: 
		return 0 #If for whatever reason there is no region information available (because the mask is empty and a previous catch didn't get it)

	return maj_axis_len

def get_pred_shape_info(pred_path: Path, 
                  save_csv: bool = True): 
    '''
    From all of the predicted masks, get the volumes associated and store in a dataframe

    Parameters 
    ----------
    pred_path: Path 
        Contains all prediction npz files 
    save_csv: bool
        Whether or not to export the dataframe as a csv. Will save to the same folder as the results 

    Returns
    ----------
    shape_df: pd.DataFrame 
        Contains the mask path and the corresponding voxel volume and pixel diameter. 
    '''

    for file in pred_path.iterdir(): 
        # Check if file is an npz, if not skip
        if not str(file).endswith('.npz'): 
            continue 
        curr_npz = np.load(file)

        #Get segmentation and spacing 
        pred_seg = curr_npz['segs'] 
        spacing = curr_npz['spacing']

        #Calculate voxel volume of whole segmentation 
        univ_vox_vol = spacing[0] * spacing[1] * spacing[2]

        vox_vol = get_vox_vol(mask_array = pred_seg, 
                              unitary_vol = univ_vox_vol)
        
        #Get RERECIST diameter 
        seg_img = sitk.GetImageFromArray(pred_seg)
        resamp_seg = resample(image = seg_img, 
					spacing = 1, 
					interpolation = 'nearest')
        resamp_seg = sitk.Cast(resamp_seg, sitk.sitkUInt8)
        resamp_seg_arr = sitk.GetArrayFromImage(resamp_seg)
        mid_slice = locate_centre_slice(resamp_seg_arr)
        diam = get_major_axis(mid_slice) 
        
        # Put into dataframe 
        cols = ['filename', 
                'predVoxVol', 
                'predDiamPix']
        
        data = [file.split["/"][-1], 
                vox_vol, 
                diam]
        
        curr_df = pd.DataFrame([data], columns = cols)
        if 'shape_df' not in locals(): 
            shape_df = curr_df 
        else: 
            shape_df = pd.concat([shape_df, curr_df], ignore_index = True).reset_index(drop = True) 
        
    if save_csv: 
        shape_df.to_csv(pred_path / 'pred_seg_shape_info.csv')
    
    return shape_df           

def merge_dataset_results(metrics_df: pd.DataFrame, 
                          sim_df: pd.DataFrame, 
                          pred_shape_df: pd.DataFrame): 
    '''
    Merge the calculated metrics, simulated information, and predicted volumes into one dataframe based on filename (assumes filenames are matching one to one)

    Parameters
    ----------
    metrics_df: pd.DataFrame
        Contains all metrics calculated from the predicted volumes 
    sim_df: pd.DataFrame 
        Contains all information pertaining to the simulation data and the RECIST categories 
    pred_shape_df: pd.DataFrame 
        Contains all shape info on the predicted mask 
    '''
    # Get common information for merging (will be sampleID_rtstructID)

    metrics_df['commonFilename'] =  metrics_df['filename'].apply(common_name_finder)
    sim_df['commonFilename'] = sim_df['Filename'].apply(common_name_finder) 
    pred_shape_df['commonFilename'] = pred_shape_df['filename'].apply(common_name_finder)

    combined_df = reduce(lambda x, y: pd.merge(x, y, on = 'commonFilename', how = 'outer'), [metrics_df, sim_df, pred_shape_df])

    # Drop any rows that contain NaNs 
    combined_df = combined_df.dropna()

    return combined_df 

def get_recist_cat(diam_change:float):
	'''
	Get the RECIST category from the diameter change. 

	Parameters
	----------
	diam_change: float 
		The diameter change of from the first to second timepoint masks 
	
	Returns 
	----------
	recist_cat: str 
		Either PD, SD, or PR
	'''
	if diam_change <= -30: 
		return 'PR'
	elif diam_change > -30 and diam_change < 20: 
		return 'SD' 
	else:
		return 'PD'
      
def get_all_recist_info(combined_df: pd.DataFrame):
    '''
    Get the predicted percentage change in voxel volume and diameter as well as the new RECIST category

    Parameters
    ----------
    combined_df: pd.DataFrame 
        Contains all information for the predictions and simulations so far

    Returns 
    ----------
    combined_df: pd.DataFrame 
        The same as combined_df but has two extra columns for the predicted RERECIST measurement, the change in volume, and information about the new RECIST category
    '''
    combined_df['PredPercChange'] = (combined_df['T1DiamPix'] - combined_df['predDiamPix']) / combined_df['predDiamPix'] * 100
    combined_df['PredVoxVolPercChange'] = (combined_df['T1VoxVol'] - combined_df['predVoxVol']) / combined_df['predVoxVol'] * 100
    combined_df['PredRECISTCat'] = combined_df['PredPercChange'].apply(get_recist_cat)

    return combined_df

def dataset_info_merging(dataset: str, 
                         disease_site: str, 
                         model_type: str, 
                         save_path: Path): 
    sim_path = dirs.PROCDATA / disease_site / dataset / 'images' / Path('sim_t1_' + dataset) / 'sim_t1_masks_summary.csv'
    metrics_path = dirs.RESULTS / disease_site / dataset / 'predictions' / model_type / 'metric_eval' / 'medsam2recist_seg_eval.csv'
    pred_path = dirs.RESULTS / disease_site / dataset / 'predictions' / model_type

    sim_df = pd.read_csv(sim_path) 
    metrics_df = pd.read_csv(metrics_path) 
    pred_df = get_pred_shape_info(pred_path)

    combined_df = merge_dataset_results(metrics_df = metrics_df, 
                                        sim_df = sim_df, 
                                        pred_shape_df = pred_df)
    
    new_combined_df = get_all_recist_info(combined_df = combined_df)

    new_combined_df.to_csv(save_path / Path('indiv_' + dataset + '_results.csv'))

def merge_all_dataset_info(all_data_path: Path): 
    '''
    Concatenate all datasets together to get one final results csv

    Parameters
    ----------
    all_data_path: Path 
        Where all of the individual csvs are located. Assumes they all have the prefix 'indiv_'.
    '''
    for file in all_data_path.iterdir():
        curr_filename = str(file).split('/')[-1] #Gets just the csv name 
        if curr_filename.startswith('indiv_'): 
            #Load csv
            curr_df = pd.read_csv(file) 

            if 'all_data_df' not in locals(): 
                all_data_df = curr_df
            else: 
                all_data_df = pd.concat([all_data_df, curr_df], ignore_index = True).reset_index(drop = True) 

    all_data_df.to_csv(all_data_path / 'combined_results.csv')

@click.command() 
@click.option('--dataset') 
@click.option('--disease_site') 
@click.option('--model_type')
@click.option('--save_path') 
def run_merging(dataset: str, 
                disease_site: str, 
                model_type: str, 
                save_path: Path): 
    dataset_info_merging(dataset = dataset, 
                         disease_site = disease_site, 
                         model_type = model_type, 
                         save_path = save_path)
    # merge_all_dataset_info(all_data_path = save_path)

if __name__ == "__main__": 
    run_merging()