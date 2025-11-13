import pandas as pd
import numpy as np
import SimpleITK as sitk
import re
import ast 
import click

from tqdm import tqdm
from joblib import Parallel, delayed
from pathlib import Path
from skimage.draw import line 
from damply import dirs

def apply_windowing(img_array: np.ndarray,
                    window_level: int, 
                    window_width: int):
    '''
    Window an image based on a window width (width of range of values to use) and a window level (where to center a window level). Otherwise known as clipping or clamping in image processing.
    
    Parameters
    ----------
    img_array: np.ndarray, 
        The image to be windowed 
    window_level: int
        Where to center the range defined in window_width
    window_width: int 
        How wide the of a range to include, centered on the level.  

    Returns 
    ----------
    windowed_img: np.ndarray
        The processed image with values clamped at the upper and lower value
    '''
    #Calculate upper and lower clamp values
    upper_val = window_level + window_width / 2 
    lower_val = window_level - window_width / 2 

    #Window image
    windowed_img = np.clip(img_array, lower_val, upper_val)

    return windowed_img 

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
    #Generate an array in the same size as the image filled with all zeros 
    recist_arr = np.zeros((img_size[0], img_size[1], img_size[2]), dtype = int)
    
    #Round the coordinate values to their nearest integers 
    coords_round = np.rint(recist_coords).astype(int)

    #Draw line using coordinates 
    rr, cc = line(coords_round[0], coords_round[1], coords_round[2], coords_round[3])

    #Put line into the correct slice in the RECIST array of all zeros 
    recist_arr[slice_number][cc, rr] = 1

    return recist_arr

def combine_bbox_match_data(matched_data: pd.DataFrame, 
                            bbox_data: pd.DataFrame): 
    '''
    Merge the bounding box data and matching annotation-image-segmentation data by the 
    annotation series instance UID and remove rows that have None values (no match found).

    Parameters
    ----------
    matched_data: pd.DataFrame
        Contains all of the data for the annotation-image-segmentation matches
    bbox_data: pd.DataFrame
        Contains all of the data for the bounding boxes from the structured reports

    Returns
    ---------
    match_bbox_data: pd.DataFrame 
        Contains all matched bounding box and annotation-image-segmentation matching data
    '''
    # Match data type for measurements 
    bbox_data['meta_MeasurementLength'] = bbox_data['meta_MeasurementLength'].astype(float)

    #Merge on the RECIST annotation series instance UID 
    merged_df = pd.merge(left = matched_data, 
                         right = bbox_data, 
                         how = 'left', 
                         left_on = ['AnnSeriesInstanceUID', 'AnnLongAxisLength'], 
                         right_on = ['meta_SeriesInstanceUID', 'meta_MeasurementLength'])
    
    #Clean dataframe to only have rows without missing values 
    match_bbox_data = merged_df.dropna()

    return match_bbox_data

def create_npzs(curr_data: pd.DataFrame, 
                save_path: Path, 
                window_width: int = None, 
                window_level: int = None):
    '''
    From the matched data and bounding box data, create .npz files for each patient and save them out to the save path specified. 

    Parameters
    ----------
    curr_data: pd.DataFrame
        Contains all matched bounding box and annotation-image-segmentation data for a given row
    save_path: Path
        Where you want the .npz files to save to 
    window_width: int 
        The width of the windowing interval. 
    window_level: int
        Where you want the windowing interval to be centered. 
    ''' 
    # Create and export npz 
    if not save_path.exists(): 
        Path(save_path).mkdir(parents = True, exist_ok = True)

    #Get all other information found in the matched bbox dataframe that is necessary for .npz file creation 
    img_file = dirs.PROCDATA / Path('Abdomen/TCIA_CPTAC-CCRCC/images') / curr_data['ImgNIFTIShortLoc']
    seg_file = dirs.PROCDATA / Path('Abdomen/TCIA_CPTAC-CCRCC/images') / curr_data['SegNIFTIShortLoc']
    bbox = np.array(ast.literal_eval(curr_data['BoundingBox']))
    slice_num = curr_data['NIFTISliceBBox']
    save_name = "_".join(str(curr_data['SegNIFTIShortLoc']).split("/")[-3:-1]) # Gets the med-imagetools unique file name from the total path. 

    #Load in nifti image and segmentation as arrays 
    curr_img = sitk.ReadImage(img_file, outputPixelType = sitk.sitkInt16)
    curr_seg = sitk.ReadImage(seg_file, outputPixelType = sitk.sitkInt8)

    curr_img_array = sitk.GetArrayFromImage(curr_img)
    curr_seg_array = sitk.GetArrayFromImage(curr_seg)

    #Window image if values were entered 
    if window_width is not None: 
        curr_img_array = apply_windowing(img_array = curr_img_array, 
                                        window_level = window_level, 
                                        window_width = window_width)
    
    # Get the spacing, direction, and origin of the image 
    spacing = curr_img.GetSpacing()
    direction = curr_img.GetDirection()
    origin = curr_img.GetOrigin()
    
    # Create the RECIST line array 
    recist_ann_arr = get_line_from_recist(recist_coords = bbox, 
                                            slice_number = slice_num, 
                                            img_size = curr_img_array.shape)

    # Create unique save path
    counter = 0
    full_save_path = save_path / Path(save_name + "_" + str(counter) + ".npz")
    while full_save_path.exists():
        counter += 1 
        new_save_name = save_path / Path(save_name + "_" + str(counter) + ".npz")
        full_save_path = new_save_name
    
    np.savez_compressed(full_save_path, 
                imgs = curr_img_array, 
                gts = curr_seg_array, 
                recist = recist_ann_arr,
                spacing = spacing,
                direction = direction, 
                origin = origin)
    
    #Free up room 
    del curr_img
    del curr_seg
    del curr_img_array 
    del curr_seg_array

@click.command()
@click.option('--dataset')
@click.option('--disease_site')
@click.option('--num_jobs')
def run_npz_create(dataset: str, 
                   disease_site: str, 
                   num_jobs: int): 
    '''
    Run npz file creation.
    '''        
    matched_df_path = dirs.PROCDATA / Path(disease_site + "/" + dataset) / Path('metadata/annotation_seg_matching/matching_ann_to_seg.csv')
    bbox_df_path = dirs.PROCDATA / Path(disease_site + "/" + dataset) / Path('metadata/bbox_gen/recist_bbox_info.csv')
    jobs = num_jobs

    matched_df = pd.read_csv(matched_df_path)
    bbox_df = pd.read_csv(bbox_df_path)

    combined_data_df = combine_bbox_match_data(matched_data = matched_df, 
                                               bbox_data = bbox_df)
    
    save_path = dirs.PROCDATA / Path(disease_site + "/" + dataset) / 'images/npz_' + dataset

    #Make sure that the save path exists 
    if not save_path.exists(): 
        Path(save_path).mkdir(parents = True, exist_ok = True)

    Parallel(n_jobs = jobs)(delayed(create_npzs)
                            (curr_data = curr_data, 
                             save_path = save_path, 
                             window_level = 40, 
                             window_width = 350)
                             for _, curr_data in tqdm(combined_data_df.iterrows(), total = combined_data_df.shape[0]))
if __name__ == '__main__': 
    run_npz_create()
    ## For debugging ##
    # matched_df_path = dirs.PROCDATA / 'Abdomen/TCIA_CPTAC-CCRCC/metadata/annotation_seg_matching/matching_ann_to_seg.csv'
    # bbox_df_path = dirs.PROCDATA / 'Abdomen/TCIA_CPTAC-CCRCC/metadata/bbox_gen/recist_bbox_info.csv'
    # jobs = 4

    # matched_df = pd.read_csv(matched_df_path)
    # bbox_df = pd.read_csv(bbox_df_path)

    # combined_data_df = combine_bbox_match_data(matched_data = matched_df, 
    #                                            bbox_data = bbox_df)
    
    # save_path = dirs.PROCDATA / 'Abdomen/TCIA_CPTAC-CCRCC/bbox_gen'

    # Parallel(n_jobs = jobs)(delayed(create_npzs)
    #                         (curr_data = curr_data, 
    #                          save_path = save_path, 
    #                          window_level = 40, 
    #                          window_width = 350)
    #                          for _, curr_data in tqdm(combined_data_df.iterrows(), total = combined_data_df.shape[0]))