import SimpleITK as sitk
import numpy as np
import math
import logging
import datetime
import pandas as pd
import random 
import click 

from pathlib import Path
from skimage.morphology import binary_erosion, binary_dilation, ball
from scipy import ndimage as ndi
from skimage.measure import regionprops, label
from damply import dirs

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
	maj_axis_len = region_info[0].axis_major_length

	return maj_axis_len

def make_post_mask(mask_array: np.ndarray,
				   diam_change_bound: float, 
				   iter_start: int = 1):
	"""
	From a given mask array, dilate or erode the 3D mask until it reaches the volume change specified.

	Parameters
	----------
	mask_array: np.ndarray,
	    Contains the original 3D binary mask for a tumour
	mask_spacing: 

	diam_change_bound: float
	    The the upper/lower limit of percentage diameter change that is requested. Keep value bounded by [-80, 100] 
	iter_start: int
	    How many iterations to start with (recommended to keep at the default 1 to prevent issues with eroding or dilating too quickly)

	Returns
	----------
	post_mask_array: np.ndarray,
	    The dilated or eroded mask
	actual_perc_diam_change: float 
		The actual percentage diameter change between the original mask and the simulated mask 
	start_diam: float
		The diameter of the original mask at the middle slice. Based on pixels, not mm
	curr_diam: float
		The diameter of the simulated mask at the middle slice (same slice as start_diam). Based on pixels, not mm
	middle_slice: int
		The NIFTI slice that the diameters are calculated on	
	"""
	# Define footprint
	footprint = ndi.generate_binary_structure(
		rank=3, connectivity=2
	)  # Generates sphere

	# Calculate starting diameter
	middle_slice = locate_centre_slice(mask_array)
	start_diam = get_major_axis(mask_array[middle_slice])
	logger.info("Starting major axis length: %s", str(start_diam))
	perc_diam_change = 0  # Initialize tracker for percentage volume change

	# Choose dilation or erosion by sign of volume change
	if diam_change_bound > 0:
		logger.info('Dilation selected for a diameter change bound of %s', str(diam_change_bound))
		first_flag = 1
		while perc_diam_change <= diam_change_bound:
			logger.info('Performing mask dilation')
			curr_mask = ndi.binary_dilation(mask_array, footprint, iterations=iter_start)  # Perform the dilation
			logger.info('Dilation complete')
			curr_diam = get_major_axis(curr_mask[middle_slice])
			logger.info("Current diameter: %s", str(curr_diam))
			perc_diam_change = (curr_diam - start_diam) / start_diam * 100
			logger.info('Current percent diameter change: %s', str(perc_diam_change))
			if (perc_diam_change > diam_change_bound):  #If percentage change is now too big
				if first_flag: #If the mask is too small to start, the dilation will not be able increase by fine amounts. Will return the current dilated mask.
					logger.debug("Too small volume to dilate precisely. Returning the original mask")
					return mask_array, 0, start_diam, start_diam, middle_slice
				else: 
					logger.info("Dilation now over upper bound. Returning the previous dilated mask")
					return prev_mask, prev_perc_diam_change, start_diam, prev_diam, middle_slice  #Take the previous mask iteration if this was not the first iteration
				
			elif perc_diam_change < diam_change_bound:
				logger.info('Diameter change stil under the dilation upper bound. Continuing to iterate')
				prev_mask = curr_mask
				prev_diam = curr_diam
				prev_perc_diam_change = perc_diam_change
				first_flag = 0
				iter_start += 1

	# Negative diameter changes imply erosion
	elif diam_change_bound < 0:
		logger.info('Erosion selected for a diameter change bound of %s', str(diam_change_bound))
		first_flag = 1
		while perc_diam_change >= diam_change_bound:
			logger.info('Performing mask erosion')
			curr_mask = ndi.binary_erosion(mask_array, footprint, iterations=iter_start)  # Perform the erosion
			logger.info('Erosion complete')
			# Check if erosion has removed all mask 
			if np.count_nonzero(curr_mask[middle_slice].astype(int)) == 0: # If there are no non-zero values in the current mask 
				logger.debug("Erosion has caused disappearance of all non-zero values in mask. Returning previous mask.")
				if first_flag: 
					return mask_array, 0, start_diam, start_diam, middle_slice
				else: 
					return prev_mask, prev_perc_diam_change, start_diam, prev_diam, middle_slice
			curr_diam = get_major_axis(curr_mask[middle_slice])
			logger.info("Current diameter: %s", str(curr_diam))
			perc_diam_change = (curr_diam - start_diam) / start_diam * 100
			logger.info('Current diameter change: %s', str(perc_diam_change))
			if perc_diam_change < diam_change_bound: 
				if first_flag: #If the mask is too small to start, the erosion will not be able decrease by fine amounts. Will return original mask 
					logger.debug("Too small volume to erode precisely. Returning original mask")
					return mask_array, 0, start_diam, start_diam, middle_slice
				else: 
					logger.info("Erosion now under lower bound. Returning the previous eroded mask")
					return prev_mask, prev_perc_diam_change, start_diam, prev_diam, middle_slice #Return the previous mask iteration if this was not the first iteration

			elif perc_diam_change > diam_change_bound:
				logger.info('Diameter change is still above the lower bound. Continuing to iterate.')
				prev_mask = curr_mask
				prev_diam = curr_diam
				prev_perc_diam_change = perc_diam_change
				first_flag = 0
				iter_start += 1


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

click.command() 
click.option('--dataset')
click.option('--disease_site')
click.option('--save_folder_name')
def run_sim_second_timepoint(dataset: str, 
							 disease_site: str, 
							 save_folder_name: Path, 
							 get_summary_csv: bool = True): 
	'''
	From ground truth segmentation NIFTI files, create a dilated or eroded version to mimic a treatment response, stable disease, or 
	progression of disease. 

	Parameters
	----------
	dataset: str
		The full name of the dataset being used (corresponds to the folder name in the common data structure) 
	disease_site: str
		The anatomical area where the disease is located (corresponds to the folder name in the common data structure)
	save_folder_name: Path
		The name of the folder to store the modified masks in. Will be stored in data/procdata/<DISEASE_SITE>/<DATASET>/images
	get_summary_csv: bool 
		Save a summary of all of the dilations/erosions done on which files, how much of a change there was in terms of diameter 
		and volume, and the simulated RECIST category. Default is true. Will be saved to save_folder_name. 
	'''
	# Get relevant paths 
	dataset_short = dataset.split("_")[-1] 
	idx_nifti_path = dirs.PROCDATA / disease_site / dataset / Path("images/mit_" + dataset_short) / Path("mit_" + dataset_short + "_index-simple.csv")
	out_path = dirs.PROCDATA / disease_site / dataset / Path("images") / save_folder_name
	
	if not out_path.exists():
		out_path.mkdir(parents=True, exist_ok = True)

	cols = ['Filename', 
		 	'PercChangeRequested', 
		 	'ActualPercChange', 
		 	'RECISTCat',
		 	'T0DiamPix', 
		 	'T1DiamPix', 
		 	'MidSlice',
		 	'T0VoxVol', 
		 	'T1VoxVol', 
		 	'VoxVolPercChange']
	
	# Go through all that were processed 
	nifti_index = pd.read_csv(idx_nifti_path) 
	seg_only_index = nifti_index[(nifti_index['Modality'] == 'RTSTRUCT') | (nifti_index['Modality'] == 'SEG')]

	for _, row in seg_only_index.iterrows(): 
		rel_filepath = row['filepath']
		full_filepath = dirs.PROCDATA / disease_site / dataset / Path("images/mit_" + dataset_short) / Path(rel_filepath)

		# Load in the NIFTI segmentation and get the relevant information for conversion 
		nifti_seg = sitk.ReadImage(full_filepath)
		nifti_mask = sitk.GetArrayFromImage(nifti_seg)

		orig_spacing = nifti_seg.GetSpacing()
		orig_direction = nifti_seg.GetDirection()
		orig_origin = nifti_seg.GetOrigin()
		unitary_vox_vol = orig_spacing[0] * orig_spacing[1] * orig_spacing[2]
		# Choose dilation or erosion at random 
		diam_change_bound = 0
		while diam_change_bound == 0: #Prevents no change 
			diam_change_bound = random.randint(-99, 100)

		t1_mask, perc_change, start_diam_pix, t1_diam_pix, mid_slice = make_post_mask(mask_array = nifti_mask, 
																				diam_change_bound = diam_change_bound)
		
		# Get t0 and t1 volumes 
		t0_vol = get_vox_vol(nifti_mask, unitary_vox_vol)
		t1_vol = get_vox_vol(t1_mask, unitary_vox_vol)
		vol_change = (t1_vol - t0_vol) / t0_vol * 100

		# Give generated mask all necessary properties before exporting to NIFTI 
		t1_img = sitk.GetImageFromArray(t1_mask.astype(int)) 

		t1_img.SetSpacing(orig_spacing)
		t1_img.SetDirection(orig_direction) 
		t1_img.SetOrigin(orig_origin) 

		img_savename = "_".join(rel_filepath.split("/")[:2])
		img_savepath = out_path / Path(img_savename + ".nii.gz")
		sitk.WriteImage(t1_img, img_savepath)

		# Get simulated RECIST category 
		recist_cat = get_recist_cat(perc_change)

		# Write intermediate data to pandas dataframe 
		if get_summary_csv: 
			curr_data = [rel_filepath, 
						diam_change_bound, 
						perc_change, 
						recist_cat,
						start_diam_pix, 
						t1_diam_pix, 
						mid_slice, 
						t0_vol, 
						t1_vol, 
						vol_change]
			
			curr_df = pd.DataFrame([curr_data], columns = cols, index = [0]) 
			if 'sim_t1_masks_info' not in locals(): 
				sim_t1_masks_info = curr_df
			else: 
				sim_t1_masks_info = pd.concat([sim_t1_masks_info, curr_df], ignore_index = True).reset_index(drop = True) 
	
	if get_summary_csv: 
		sim_t1_masks_info.to_csv(out_path / 'sim_t1_masks_summary.csv')

if __name__ == '__main__': 
	# Get current time for logging 
	curr_datetime = datetime.datetime.now()
	str_datetime = curr_datetime.strftime('%Y-%m-%d %H:%M:%S')
	#Configure logger
	logger = logging.getLogger(__name__)
	logging.basicConfig(filename = dirs.LOGS / Path("sim_second_timepoint_" + str_datetime + ".log"), encoding='utf-8', level=logging.DEBUG)

	run_sim_second_timepoint(dataset = 'TCIA_CPTAC-CCRCC', 
						  disease_site = 'Abdomen', 
						  save_folder_name = 't1_segs_test')
							


	### PREVIOUS TESTING ###
	# mask_img = sitk.ReadImage(test_mask)
	# mask_spacing = mask_img.GetSpacing()
	# print(mask_spacing)
	# vox_unit_vol = mask_spacing[0] * mask_spacing[1] * mask_spacing[2]

	# mask_arr = sitk.GetArrayFromImage(mask_img)

	# ## Testing dilation
	# perc_grow = -40

	# dilated_mask = make_post_mask(mask_array=mask_arr, 
	# 						diam_change_bound=perc_grow)

	# start_vol = get_vox_vol(mask_arr, vox_unit_vol)
	# print(start_vol)
	# dilated_vol = get_vox_vol(dilated_mask, vox_unit_vol)
	# print(dilated_vol)

	# actual_perc_change = (dilated_vol - start_vol) / start_vol * 100

	# print('Percentage Volume Change: ', actual_perc_change)

	# # Get image from new mask 
	# # Get center slice 
	# center = locate_centre_slice(mask_arr)
	# dilated_slice = dilated_mask[center]
	# non_dil_slice = mask_arr[center]

	# region_dil = regionprops(dilated_slice.astype(int)) 
	# major_axis = region_dil[0].axis_major_length
	# print("Major Axis for dilated mask: ", major_axis)

	# region_nondil = regionprops(non_dil_slice.astype(int))
	# major_axis_nondil = region_nondil[0].axis_major_length
	# print("Major Axis for non-dilated mask: ", major_axis_nondil)

	# actual_perc_diam_grow = (major_axis - major_axis_nondil) / major_axis_nondil * 100

	# print("Major Axis percentage diameter change: ", actual_perc_diam_grow)

	# dilated_img = sitk.GetImageFromArray(dilated_mask.astype(int))
	# # Set direction, spacing, and orientation to be the same as the input array 
	# dilated_img.SetSpacing(mask_img.GetSpacing())
	# dilated_img.SetDirection(mask_img.GetDirection())
	# dilated_img.SetOrigin(mask_img.GetOrigin())
	# sitk.WriteImage(dilated_img, 'test_erosion_mask_sphere40_C3L-00792_0001_RTSTRUCT_617561.4.nii.gz')