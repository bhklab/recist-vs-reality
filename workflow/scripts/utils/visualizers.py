import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path
from skimage.measure import label, regionprops


def locate_centre_slice(mask_3d: np.ndarray):
    """
    Locates the center slice of a 3D mask.

    Args:
    - mask_3d (ndarray): A 3D binary mask.

    Returns:
    - (int): The index of the center slice.
    """
    
    # find the slice in the center
    lesion_labels = label(mask_3d)
    centre_slc = regionprops(lesion_labels)[0].centroid[2]

    # number of voxels per slice
    vox_per_slc = np.array([np.sum(mask_3d[:,:,slc]) for slc in range(mask_3d.shape[2])])
    max_vox_slc = np.where(vox_per_slc==np.max(vox_per_slc))[0]

    return int(np.floor((centre_slc+max_vox_slc)/2))


def get_image_range(window_level: float | int,
                    window_width: float | int):
    
    upper_grey_level = np.ceil(window_level + (window_width / 2))
    lower_grey_level = np.floor(window_level - (window_width / 2))

    return lower_grey_level, upper_grey_level


def visualize_2D_segmentations(image: np.ndarray, 
                               ground_truth_mask: np.ndarray, 
                               predicted_mask: np.ndarray | None = None, 
                               slice_index: int | None = None,
                               alpha: float = 0.3,
                               window_level: float | int | None = None,
                               window_width: float | int | None = None,
                               save_path: str | Path | None = None,
                               sample_id: str | Path | None = None, 
                               dice: float | None = None, 
                               hd_95: float | None = None, 
                               apl: float | None = None, 
                               input_box=None):
    """
    Visualize 2D segmentations.
    
    Args:
    - image (np.ndarray): The image to visualize.
    - ground_truth_mask (np.ndarray): The ground truth segmentation mask.
    - predicted_mask (np.ndarray): The predicted segmentation mask.
    - alpha (float): Value to set transparency for mask overlays
    - slice_index (int): The index of the slice to plot.
    - save_path (str, Path): Path to save the figure to, including file name and extension.
    - dice (float): DICE score to include on the plot as a subtitle.
    - hd95 (float): Hausdorff Distance 95th Percentile score to include on the plot as a subtitle.
    - apl (float): Added Path Length score to include on the plot as a subtitle.
    - input_box (): Bounding box used as input for segmentation model, to plot on figure.

    Returns:
    - plt.Figure: Indicated slice (or center slice of segmentation) of image with ground truth and predicted masks overlaid
    """
    # If no slice index provided, get the centre slice of the ROI in the ground truth mask
    if not slice_index:
        slice_index = locate_centre_slice(ground_truth_mask)

    fig, ax = plt.subplots(figsize=(8, 8))

    if window_width is None or window_level is None:
        lower_grey_level = image.min()
        upper_grey_level = image.max()
    else:
        lower_grey_level, upper_grey_level = get_image_range(window_level, window_width)

    # Plot slice of CT
    ax.imshow(
        image[:, :, slice_index], 
        cmap=plt.cm.Greys_r, 
        vmin=lower_grey_level,
        vmax=upper_grey_level
    )

    # Make mask of original ROI in green
    ground_truth_mask_mask = np.ma.masked_where(ground_truth_mask == 0, ground_truth_mask)
    ax.imshow(
        ground_truth_mask_mask[:, :, slice_index],
        cmap=plt.cm.RdYlGn,
        vmin=ground_truth_mask_mask.min(),
        vmax=ground_truth_mask_mask.max(),
        alpha=alpha,
    )

    # Make mask of predicted segmentation image
    if predicted_mask is not None:
        # If segmentation prediction was 3D, pick out the indicated slice prior to plotting
        if predicted_mask.ndim == 3:
            predicted_mask = predicted_mask[:, :, slice_index]
        
        predicted_mask_mask = np.ma.masked_where(predicted_mask == 0, predicted_mask)
        ax.imshow(
            predicted_mask_mask[:, :],
            cmap=plt.cm.cool,
            vmin=predicted_mask_mask.min(),
            vmax=predicted_mask_mask.max(),
            alpha=alpha,
        )

    # If bounding box provided, plot it
    if input_box is not None:
        x_min, y_min, x_max, y_max = input_box[0] 

        # Input box is 1024 px, Segmentations are varied numbers of px
        scale_x = ground_truth_mask.shape[1] / 1024
        scale_y = ground_truth_mask.shape[0] / 1024

        x_min *= scale_x
        x_max *= scale_x
        y_min *= scale_y
        y_max *= scale_y

        rect = patches.Rectangle(
            (x_min, y_min), 
            x_max - x_min,  
            y_max - y_min, 
            linewidth=1,
            edgecolor='red',
            facecolor='none'
        )
        ax.add_patch(rect)

    ax.axis("off")
    
    if sample_id:
        plt.suptitle(sample_id)

    # Add metrics
    caption = ""
    if dice is not None:
        caption += f'DSC: {dice:.3f}'
    if hd_95 is not None:
        caption += f', 95HD: {hd_95:.2f}'
    if apl is not None:
        caption += f', APL: {apl:.2f}'

    if caption:
        ax.set_title(caption)
    # fig.text(0.5, 0.0, caption, ha='center', va='center')

    if save_path:
        # Save figure
        plt.savefig(save_path, bbox_inches='tight', pad_inches=0.1)
        plt.close()
    
    return plt
