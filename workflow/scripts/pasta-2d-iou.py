from damply import dirs
import pandas as pd
import SimpleITK as sitk
import numpy as np
from skimage.measure import label
from imgtools.transforms.functional import resample

import logging

print(dirs)
dirs.LOGS.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
	level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
	filename=dirs.LOGS / "pasta_2d_iou.log"
)

logger = logging.getLogger(__name__)

def load_img(path,
             desired_spacing=[1.0,1.0,1.0]
             ) -> tuple[np.array, int]:

    logger.info(f"Loading {path}...")

    img = sitk.ReadImage(path)
    img_spacing = list(img.GetSpacing())

    if img_spacing != desired_spacing:
        logger.info(f"Image was not isotropic, respacing to 1x1x1.")
        img = resample(img, spacing = desired_spacing, interpolation='nearest')
        
        img_arr = sitk.GetArrayFromImage(img)
        # Something about the resampling flips the image in the sagittal and coronal planes
        logger.info("Flipping mask back in sagittal and coronal planes")
        img_arr = np.flip(img_arr, axis=[1,2])

    else:
        img_arr = sitk.GetArrayFromImage(img)

    if np.count_nonzero(img_arr) > 0: 
        conn_comps, num_comps = label(img_arr, return_num=True)
        return conn_comps, num_comps
    
    else:
        logger.info(f"Mask file at {path} is empty. No predictions made.")
        return img_arr, 0



def get_max_slice(np_mask:np.array) -> tuple[np.array, int]:
    # Sum the mask in the x and y axes to find the axial slice with the largest tumour area
    axial_sum = np.sum(np_mask, axis=(1,2))
    # Get the index of the axial slice with the largest tumour area
    max_axial_index = np.argmax(axial_sum)
    # Select out the slice with the largest index
    max_area_slice = np_mask[max_axial_index]

    return max_area_slice, max_axial_index


def bbox_from_seg(seg_2d: np.ndarray) -> list[np.int64]: 
    '''  
    Get a bounding box from a given segmentation (binary array). Assumes that only one tumour is within this segmentation.

    Parameters
    ----------
    seg2d: np.ndarray
        The segmentation slice to get a bounding box from 

    Returns
    ----------
    bbox: list 
        A list of coordinates for the bounding box in the order of [xmin, ymin, xmax, ymax]
    '''
    rows = np.any(seg_2d, axis=1)
    cols = np.any(seg_2d, axis=0)
    ymin, ymax = np.where(rows)[0][[0, -1]]
    xmin, xmax = np.where(cols)[0][[0, -1]]

    bbox = [xmin, ymin, xmax, ymax]
    
    return bbox


def proc_component(comps:np.ndarray,
                   comp_idx:int,
                   thresh_inc:int = 10
                   ) -> tuple[np.array, int] | tuple[int, int]:
    # Set all other components to 0s
    curr_comp = np.ma.masked_where(comps != comp_idx, comps).filled(0)

    # check if size of component is large enough
    if np.count_nonzero(curr_comp) < thresh_inc:
            return -1, -1
      
    curr_comp_bin = np.clip(curr_comp, 0, 1)
    curr_comp_bin = curr_comp_bin.astype(int)

    _max_axial_slice, max_axial_index = get_max_slice(curr_comp_bin)

    return curr_comp_bin, max_axial_index


def calc_2d_IoU(bbox_gt: np.ndarray,
                bbox_pred: np.ndarray
                ) -> np.float64:
    '''  
    Calculate the 2D intersection-over-union (IoU) between two bounding boxes. 

    Parameters
    ----------
    bbox_gt: np.ndarray
        Ground truth 2D bounding box in the order [xmin, ymin, xmax, ymax]
    bbox_pred: np.ndarray
        Predicted 2D bounding box in the order [xmin, ymin, xmax, ymax]

    Returns 
    ----------
    iou: float 
        The IoU of the two bounding boxes
    '''
    # Get coordinates and area of intersecting rectangle
    xmin_inter = max(bbox_gt[0], bbox_pred[0])
    ymin_inter = max(bbox_gt[1], bbox_pred[1])
    xmax_inter = min(bbox_gt[2], bbox_pred[2])
    ymax_inter = min(bbox_gt[3], bbox_pred[3])

    inter_width = max(0, xmax_inter-xmin_inter)
    inter_height = max(0, ymax_inter-ymin_inter)
    inter_area = inter_width * inter_height 

    # Get area of the each bounding box 
    gt_area = (bbox_gt[2]-bbox_gt[0]) * (bbox_gt[3]-bbox_gt[1])
    pred_area = (bbox_pred[2]-bbox_pred[0]) * (bbox_pred[3]-bbox_pred[1])

    # Calculate union area 
    union_area = gt_area + pred_area - inter_area 

    # Calculate IoU 
    if union_area == 0: 
        return 0.0 # This should never happen, but implementing just in case
    
    iou = inter_area / union_area 

    return iou


def match_components(sample_id,
                     gt_path,
                     pred_path,
                     thresh_inc:int = 10):
    
    gt_conn_comps, gt_num_comps = load_img(gt_path)
    pred_conn_comps, pred_num_comps = load_img(pred_path)
    match_dict = {}

    if pred_num_comps == 0:
        logger.info("Empty prediction mask. No calculations done.")
        return match_dict

    logger.info("Masks loaded")
    
    for gt_idx in range(1, gt_num_comps+1):
        # Getting single ground truth ROI (component)
        curr_gt_comp, gt_max_axial_index = proc_component(gt_conn_comps, gt_idx, thresh_inc)

        # Check if ground truth ROI is too small for calculation
        if isinstance(curr_gt_comp, int) and curr_gt_comp -1:
            logger.info(f"{sample_id} ground truth ROI mask {gt_idx} is less than {thresh_inc} voxels. Skipping calculation.")
            continue

        # Check if there are any overlaps with the individual predicted segmentations 
        for pred_idx in range(1,pred_num_comps+1): 
            curr_pred_comp, _pred_max_axial_index = proc_component(pred_conn_comps, pred_idx, thresh_inc)

            if isinstance(curr_pred_comp, int) and curr_pred_comp -1:
                logger.info(f"{sample_id} predicted ROI mask {pred_idx} is less than {thresh_inc} voxels. Skipping calculation.")
                continue

            curr_intersect = curr_gt_comp & curr_pred_comp
            if np.count_nonzero(curr_intersect) > 0: 
                # Indicates overlap between the two segmentations, so these have matched. Save as a pair for further analysis
                logger.info(f"Overlapping ground truth and predicted mask found for {sample_id}")

                logger.info("Getting bounding boxes for ground truth and predicte mask...")
                gt_bbox = bbox_from_seg(curr_gt_comp[gt_max_axial_index])
                pred_bbox = bbox_from_seg(curr_pred_comp[gt_max_axial_index])

                logger.info("Calculating 2D IoU...")
                match_2D_iou = calc_2d_IoU(gt_bbox, pred_bbox)

                match_id = f"{sample_id.removesuffix(".nii.gz")}_gt{gt_idx}_pred{pred_idx}"
                match_dict[match_id] = {
                    "sample_id": sample_id,
                    "gt_bbox": list(map(int, gt_bbox)),
                    "pred_bbox": list(map(int, pred_bbox)),
                    "slice_idx_from_gt": int(gt_max_axial_index),
                    "slice_2D_iou": float(match_2D_iou) 
                }

    return match_dict


if __name__ == "__main__":
    print(dirs.LOGS)
    dirs.LOGS.mkdir(parents=True, exist_ok=True)

    dataset = "PASTA"
    gt_folder = dirs.RAWDATA / dataset / "gt"
    pred_folder = dirs.PROCDATA / dataset / "pred"
    results_folder = dirs.RESULTS / dataset 
    results_folder.mkdir(parents=True, exist_ok=True)

    gt_files = list(gt_folder.glob("*.nii.gz"))
    pred_files = list(pred_folder.glob("*.nii.gz"))

    # Get set of samples for each
    gt_samples = {file.name for file in gt_files}
    pred_samples = {file.name for file in pred_files}

    samples_to_proc = pred_samples & gt_samples

    results = {}
    for sample_id in samples_to_proc:
        sample_results = match_components(sample_id = sample_id,
                                   gt_path = gt_folder / sample_id,
                                   pred_path = pred_folder / sample_id)
        
        results.update(sample_results)
        
    results_df = pd.DataFrame().from_dict(results, orient='index')


    results_df.to_csv(results_folder / f"{dataset}_2D_IoU_results.csv", index=True, index_label='match_id')
    logger.info("Script completed.")