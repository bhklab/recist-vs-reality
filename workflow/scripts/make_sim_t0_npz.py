import numpy as np
import matplotlib.pyplot as plt
import random 
import click 

from damply import dirs
from pathlib import Path 
from skimage.draw import line

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

def elongate_shrink_recist(coords, factor=1.0):
    """
    Elongate or shrink a RECIST line by a given factor.
    
    Args:
        coords: numpy array of shape (N, 2) with line coordinates
        factor: float, factor to elongate (>1) or shrink (<1) the line
                e.g., 1.2 = 20% longer, 0.8 = 20% shorter
    
    Returns:
        numpy array of shape (N, 2) with adjusted coordinates
    """
    coords = np.array(coords)
    
    # Get endpoints
    start_point = coords[0]
    end_point = coords[-1]
    
    # Calculate center point
    center = (start_point + end_point) / 2.0
    
    # Calculate direction vector
    direction = end_point - start_point
    length = np.linalg.norm(direction)
    
    if length == 0:
        return coords  # Degenerate case
    
    # Normalize direction
    direction_unit = direction / length
    
    # Calculate new endpoints
    new_half_length = (length * factor) / 2.0
    new_start = center - direction_unit * new_half_length
    new_end = center + direction_unit * new_half_length
    
    # Generate new coordinates along the line
    # Interpolate between new_start and new_end
    num_points = len(coords)
    t_values = np.linspace(0, 1, num_points)
    new_coords = new_start[None, :] + t_values[:, None] * (new_end - new_start)[None, :]
    
    # Round to nearest integer (pixel coordinates)
    new_coords = np.round(new_coords).astype(int)
    
    return new_coords

def adjust_recist_lines(mask, image, factor=1.0, show_plot=True):
    """
    Adjust RECIST lines on a 2D slice and visualize original vs adjusted.
    
    Args:
        mask: 3D numpy array with RECIST line mask (values == 1 indicate RECIST line)
        image: 3D numpy array with image data
        factor: float, factor to elongate (>1) or shrink (<1) the line
        show_plot: bool, whether to display the visualization
    
    Returns:
        tuple: (original_coords, adjusted_coords, slice_idx)
    """
    # Find the slice containing the RECIST line
    slice_idx = np.argwhere(mask == 1)[0][0]
    mask_2D = mask[slice_idx]
    image_2D = image[slice_idx]

    print(f"Mask 2D shape: {mask_2D.shape}")
    print(f"Image shape: {image.shape}")
    print(f"Slice index: {slice_idx}")
    
    # Extract RECIST coordinates from the mask
    recist_coords = np.argwhere(mask_2D == 1)

    if len(recist_coords) == 0:
        raise ValueError("No RECIST coordinates found in mask!")

    
    # Sort coordinates to get a proper line order (by distance from first point)
    sorted_coords = [recist_coords[0]]
    remaining = list(recist_coords[1:])
    
    # Greedily find nearest neighbors to form a line
    while remaining:
        last_point = sorted_coords[-1]
        distances = [np.linalg.norm(p - last_point) for p in remaining]
        nearest_idx = np.argmin(distances)
        sorted_coords.append(remaining.pop(nearest_idx))
    
    original_coords = np.array(sorted_coords)
    
    # Adjust the RECIST line
    adjusted_coords = elongate_shrink_recist(original_coords, factor=factor)
    
    print(f"\nOriginal RECIST line:")
    print(f"  Start: {original_coords[0]}, End: {original_coords[-1]}")
    print(f"  Length: {np.linalg.norm(original_coords[-1] - original_coords[0]):.2f} pixels")
    print(f"  Number of points: {len(original_coords)}")
    
    print(f"\nAdjusted RECIST line (factor={factor}):")
    print(f"  Start: {adjusted_coords[0]}, End: {adjusted_coords[-1]}")
    print(f"  Length: {np.linalg.norm(adjusted_coords[-1] - adjusted_coords[0]):.2f} pixels")
    print(f"  Number of points: {len(adjusted_coords)}")
    
    # Create visualization
    if show_plot:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        
        # Display the image slice
        ax.imshow(image_2D, cmap='gray', alpha=0.7)
        
        # Plot original RECIST line in red
        ax.plot(original_coords[:, 1], original_coords[:, 0], 
                'r-', linewidth=2, label='Original RECIST', alpha=0.8)
        ax.scatter(original_coords[0, 1], original_coords[0, 0], 
                    c='red', s=100, marker='o', label='Original Start', zorder=5)
        ax.scatter(original_coords[-1, 1], original_coords[-1, 0], 
                    c='darkred', s=100, marker='s', label='Original End', zorder=5)
        
        # Plot adjusted RECIST line in green
        ax.plot(adjusted_coords[:, 1], adjusted_coords[:, 0], 
                'g-', linewidth=2, label=f'Adjusted RECIST (factor={factor})', alpha=0.8)
        ax.scatter(adjusted_coords[0, 1], adjusted_coords[0, 0], 
                    c='lime', s=100, marker='o', label='Adjusted Start', zorder=5)
        ax.scatter(adjusted_coords[-1, 1], adjusted_coords[-1, 0], 
                    c='darkgreen', s=100, marker='s', label='Adjusted End', zorder=5)
        
        ax.set_title(f'RECIST Line Adjustment - Slice {slice_idx}\n'
                    f'Original Length: {np.linalg.norm(original_coords[-1] - original_coords[0]):.1f}px, '
                    f'Adjusted Length: {np.linalg.norm(adjusted_coords[-1] - adjusted_coords[0]):.1f}px',
                    fontsize=12)
        ax.legend(loc='upper right')
        ax.set_xlabel('X (columns)')
        ax.set_ylabel('Y (rows)')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    return original_coords, adjusted_coords, slice_idx

@click.command()
@click.option('--dataset')
@click.option('--disease_site')
@click.option('--npz_folder')
@click.option('--out_file_name')
def run_adjust_recist_line_npz(dataset: str, 
                               disease_site: str, 
                               npz_folder: str,
                               out_file_name: str): 
    npz_path = dirs.PROCDATA / disease_site / dataset / 'images' / npz_folder 

    for npz in npz_path.iterdir():
        if not str(npz).endswith('.npz'): 
            continue 
        curr_npz = np.load(npz)

        for i in range(3):
            #Choose random scaling factor (between 1.20 and 0.8)
            scaler = random.uniform(0.8, 1.2)

            _, adjusted_coords, slice_idx = adjust_recist_lines(mask = curr_npz['recist'], 
                                                        image = curr_npz['imgs'], 
                                                        factor = scaler, 
                                                        show_plot = False)
            
            # Get adjusted line back into full array form
            new_recist_coords = [adjusted_coords[0][0], adjusted_coords[0][1], adjusted_coords[-1][0], adjusted_coords[-1][1]]
            new_recist_line = get_line_from_recist(new_recist_coords,                                                slice_idx, 
                                                curr_npz['imgs'].shape)
            
            if np.count_nonzero(new_recist_line) < 5: 
                continue 
                
            out_path = dirs.PROCDATA / disease_site / dataset / 'images' / out_file_name
            if not out_path.exists(): 
                out_path.mkdir(parents = True, exist_ok = True) 
            
            file_name = str(npz).split("/")[-1] + "_" + str(i) + ".npz"
            full_out_path = out_path / file_name
            
            np.savez_compressed(full_out_path, 
                                imgs = curr_npz['imgs'], 
                                gts = curr_npz['gts'], 
                                recist = new_recist_line, 
                                spacing = curr_npz['spacing'], 
                                direction = curr_npz['direction'], 
                                origin = curr_npz['origin'])

if __name__ == "__main__":
    run_adjust_recist_line_npz()

    # data = np.load("/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/data/procdata/Abdomen/TCIA_CPTAC-CCRCC/images/npz_TCIA_CPTAC-CCRCC/C3N-00194_0035_RTSTRUCT_210218.4_0.npz")
    
    # original_coords, adjusted_coords, slice_idx = adjust_recist_lines(
    #     data["recist"], 
    #     data["imgs"], 
    #     factor=0.5,  
    #     show_plot=False
    # )

    # print("stop here")