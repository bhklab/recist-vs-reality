import pandas as pd 
import numpy as np 
import random 
import click

from pathlib import Path
from skimage.measure import regionprops, label
from damply import dirs

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
	maj_axis_len = region_info[0].axis_major_length

	return maj_axis_len

def get_sim_length(mask, 
                   perc_change: int):
    
    center_slice = locate_centre_slice(mask.astype(int))
    long_axis_orig = get_major_axis(mask[center_slice].astype(int)) 
    long_axis_new = long_axis_orig * (1 + perc_change/100)

    return long_axis_new

def get_recist_cat(diam_change: float):
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

@click.command() 
@click.option('--dataset')
@click.option('--disease_site')
@click.option('--npz_folder')
def gen_sim_RECIST(dataset: str, 
                   disease_site: str, 
                   npz_folder: str):
    '''
    testing 

    parameters
    ----------
    dataset
    disease_site
    npz_folder
    '''
    
    cols = ['filename',
            'orig_diam_pix', 
            'new_diam_pix', 
            'recist_cat',
            'perc_change',
            'orig_vox_vol']
    
    npz_path = dirs.PROCDATA / disease_site / dataset / 'images' / npz_folder 
    for npz in npz_path.iterdir(): 
        if not str(npz).endswith('.npz'): 
            continue 
        
        curr_npz = np.load(npz) 
        spacing = curr_npz['spacing']
        unit_vox_vol = spacing[0] * spacing[1] * spacing[2]
        orig_volume = get_vox_vol(curr_npz['gts'], 
                                  unit_vox_vol)
        mid_slice = locate_centre_slice(curr_npz['gts'].astype(int))
        long_axis_orig = get_major_axis(curr_npz['gts'][mid_slice].astype(int)) 
        
        #Generate a random response 
        perc_diam_change= random.randint(-99, 100)
        recist_cat = get_recist_cat(perc_diam_change) 

        new_long_ax = get_sim_length(curr_npz['gts'], 
                                     perc_diam_change)
        
        filename = str(npz).split("/")[-1]
        
        curr_data = [filename, 
                     long_axis_orig, 
                     new_long_ax, 
                     recist_cat, 
                     perc_diam_change,
                     orig_volume]

        curr_df = pd.DataFrame([curr_data], columns = cols, index = [0]) 
        if 'sim_t1_masks_info' not in locals(): 
            sim_t1_masks_info = curr_df
        else: 
            sim_t1_masks_info = pd.concat([sim_t1_masks_info, curr_df], ignore_index = True).reset_index(drop = True)
        
    sim_t1_masks_info.to_csv(npz_path / 'sim_recist_measures.csv') 

if __name__ == "__main__":
    gen_sim_RECIST()
