import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def visualize_2D_segmentations(image, original_seg, centre_slc_ind, path, dice, hd_95, apl, medsam_seg=None, input_box=None):
    """
    Visualize 2D segmentations.
    
    Args:
    - image (np.ndarray): The image to visualize.
    - original_seg (np.ndarray): The original segmentation.
    - medsam_seg (np.ndarray): The MedSam segmentation.
    - centre_slc_ind (int): The index of the central slice to use.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot slice of CT
    ax.imshow(
        image[:, :, centre_slc_ind], 
        cmap=plt.cm.Greys_r, 
        vmin=image.min(), 
        vmax=image.max()
    )

    # Make mask of original ROI 
    original_seg_mask = np.ma.masked_where(original_seg == 0, original_seg)
    ax.imshow(
        original_seg_mask[:, :, centre_slc_ind],
        cmap=plt.cm.PuOr,
        vmin=original_seg_mask.min(),
        vmax=original_seg_mask.max(),
        alpha=0.3,
    )

    # Make mask of MedSAM ROI 
    if medsam_seg is not None:
        medsam_seg_mask = np.ma.masked_where(medsam_seg == 0, medsam_seg)
        ax.imshow(
            medsam_seg_mask[:, :],
            cmap=plt.cm.brg,
            vmin=medsam_seg_mask.min(),
            vmax=medsam_seg_mask.max(),
            alpha=0.3,
        )

    # If bounding box provided, plot it
    if input_box is not None:
        x_min, y_min, x_max, y_max = input_box[0] 

        # Input box is 1024 px, Segmentations are varied numbers of px
        scale_x = original_seg.shape[1] / 1024
        scale_y = original_seg.shape[0] / 1024

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

    # Add metrics
    caption = f'DSC: {dice:.3f}, 95HD: {hd_95:.2f}, APL: {apl:.2f}'
    fig.text(0.5, 0.0, caption, ha='center', va='center', fontsize=20)

    # Save figure
    plt.savefig(path, bbox_inches='tight', pad_inches=0.1)
    plt.close()
