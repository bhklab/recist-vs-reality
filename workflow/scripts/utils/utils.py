import numpy as np
import pandas as pd 

from evaluate import Evaluator 

def find_first_last_slice(mask): 
    '''
    Based on a 3D mask array, get the first and last slice within the array that has masked values. 

    Parameters
    ----------
    mask: 
        3D mask array 
    
    Returns 
    ----------
    first_slice: int 
        The index where the first slice of the mask is 
    last_slice: int 
        The index where the last slice of the mask is 
    '''
    axes = tuple([i for i in range(mask.ndim) if i != 0])

    slices = mask.any(axis = axes) 

    nonzero_indices = np.where(slices)[0]

    first_slice = np.amin(nonzero_indices)
    last_slice = np.amax(nonzero_indices)

    return first_slice, last_slice

def list_nonzero_seg_slices(seg): 
    nonzero_slices = []
    for slice_idx in range(seg.shape[0]): 
        if np.count_nonzero(seg[slice_idx]) > 0: 
            nonzero_slices.append(slice_idx)
    return nonzero_slices

def calc_metrics(pred_mask: np.ndarray, 
                gt_mask: np.ndarray, 
                spacing: np.ndarray,
                filename: str, 
                text_prompt: str = None) -> pd.DataFrame:
        '''
        Calculate performance metrics based on the predicted and ground truth masks and save into a dataframe. 

        Parameters
        ----------
        pred_mask: np.ndarray
            The mask that was predicted by the model. 
        gt_mask: np.ndarray
            The ground truth segmentation array. 
        spacing: np.ndarray
            The spacing associated with the ground truth mask. 
        filename: str 
            The name of the npz file that is being evaluated.
        text_prompt: str = None
            The type of text prompt used on the image (if applicable)
        
        Returns
        ----------
        metric_df: pd.DataFrame
            Contains the evaluation performance. 
        '''
        #Initialize evaluator 
        metric_eval = Evaluator() 
        metric_dict = metric_eval(preds = pred_mask, 
                                targets = gt_mask, 
                                spacing = spacing 
                                )
        
        metric_df = pd.DataFrame(metric_dict, index = [0])

        #Add columns for the range of segmentation values (both ground truth and predicted)
        first_gts, last_gts = find_first_last_slice(gt_mask)
        try: #If no mask was predicted, this will throw an error
            first_pred, last_pred = find_first_last_slice(pred_mask) 
        except ValueError: 
            print(f"Empty predicted segmentation for file: {filename}.")
            first_pred = 0
            last_pred = 0

        # Get the list of all slices that have segmentation in them for each mask 
        mask_pred_list = list_nonzero_seg_slices(pred_mask)
        gt_list = list_nonzero_seg_slices(gt_mask) 

        metric_df['GTSliceList'] = [gt_list]
        metric_df['PredSliceList'] = [mask_pred_list]

        gts_range = [first_gts, last_gts] 
        pred_range = [first_pred, last_pred] 

        metric_df['GTSliceRange'] = [gts_range]
        metric_df['PredSliceRange'] = [pred_range]
        metric_df['filename'] = filename # To ensure we can map the results back to the segmentations 

        # Get slice interval IoU 
        metric_df['SliceIoU'] = len(list(set.intersection(set(gt_list), set(mask_pred_list)))) / len(list(set.union(set(gt_list), set(mask_pred_list))))
        metric_df['MaskUniqueSlice'] = [list(set(mask_pred_list) - set(gt_list))] # Only slices that are in predicted mask and are not in ground truth mask
        metric_df['GTUniqueSlice'] = [list(set(gt_list) - set(mask_pred_list))] # Opposite of the line above
        
        if text_prompt is not None: 
             metric_df['TextPrompt'] = text_prompt
             
        return metric_df