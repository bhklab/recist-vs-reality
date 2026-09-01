import SimpleITK as sitk
import numpy as np
from pathlib import Path

def load_nifti(path_to_nifti_file: Path | str, 
               convert_to_np_array:bool = True):
    """
    Load a NIFTI file.

    Args:
    - path_to_nifti_file (str): The path to the NIFTI file.
    - convert_to_np_array (bool): Whether to convert the NIFTI file to a NumPy array. 
                                  Will rearrange the axes so slices are the z-dimension to match sitk.Image arrangement.

    Returns:
    - ct_nifti (sitk.Image or np.ndarray): The loaded NIFTI file.
    """
    # load the ct using nibabel
    ct_nifti = sitk.ReadImage(path_to_nifti_file)

    if convert_to_np_array:
        ct_nifti = sitk.GetArrayFromImage(ct_nifti)
        # fix the axis order so the slices are in the last dimension
        ct_nifti = np.moveaxis(ct_nifti, 0, -1)
        
    return ct_nifti