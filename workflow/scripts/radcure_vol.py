### Get radcure data volume

import pandas as pd 
import numpy as np 
import click
import SimpleITK as sitk 

from joblib import Parallel, delayed
from tqdm import tqdm
from damply import dirs 
from pathlib import Path 
from functools import reduce 
from skimage.measure import regionprops, label
from imgtools.transforms.functional import resample 

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

def get_pred_shape_info(file: Path): 
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
    # Check if file is an npz, if not skip
    if not str(file).endswith('.npz'): 
        return pd.DataFrame() 
    curr_npz = np.load(file)
    print("Getting volume for file: ", file)
    #Get segmentation and spacing 
    pred_seg = curr_npz['gts'] 
    spacing = curr_npz['spacing']

    #Calculate voxel volume of whole segmentation 
    univ_vox_vol = spacing[0] * spacing[1] * spacing[2]
    
    #Get RERECIST diameter 
    seg_img = sitk.GetImageFromArray(pred_seg)
    seg_img.SetSpacing(spacing)
    resamp_seg = resample(image = seg_img, 
                spacing = 1, 
                interpolation = 'nearest')
    resamp_seg = sitk.Cast(resamp_seg, sitk.sitkUInt8)
    resamp_seg_arr = sitk.GetArrayFromImage(resamp_seg)

    vox_vol = get_vox_vol(mask_array = resamp_seg_arr, 
                            unitary_vol = univ_vox_vol)
    
    # Put into dataframe 
    cols = ['filename', 
            'T0VoxVol']
    
    data = [str(file).split("/")[-1], 
            vox_vol]
    
    curr_df = pd.DataFrame([data], columns = cols)

    return curr_df 

@click.command()
@click.option("--npz_t0_path")
def run_vols(npz_t0_path: Path):
    vol_info = Parallel(n_jobs = 16)(delayed(get_pred_shape_info)
                                     (file = file)
                                     for file in tqdm(npz_t0_path.iterdir(), total = 2996))
    
    for info in vol_info: 
        if 'all_sim_info_df' not in locals(): 
            all_sim_info_df = info
        else: 
            all_sim_info_df = pd.concat([all_sim_info_df, info])

    all_sim_info_df.to_csv(dirs.RESULTS)
                    
if __name__ == "__main__": 
    run_vols()
    

    

