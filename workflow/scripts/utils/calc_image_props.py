from pathlib import Path
import SimpleITK as sitk
from skimage.measure import regionprops
import pandas as pd
import tqdm
import click
from joblib import Parallel, delayed


def ShortAxisDiameter(image: sitk.Image):
    """
    Calculate the short axis length of a binary mask image and save it to a text file.

    Parameters:
    - image_path: Path to the input binary mask image (e.g., NIfTI file).
    - output_path: Path to the output text file where the short axis length will be saved.
    """
    # Convert image to numpy array
    image_array = sitk.GetArrayFromImage(image)

    # Calculate region properties
    props = regionprops(image_array.astype(int))

    if not props:
        raise ValueError("No regions found in the binary mask.")

    # Assuming we are interested in the largest region
    largest_region = max(props, key=lambda x: x.area)

    # Get the lengths of the major and minor axes
    # major_axis_length = largest_region.major_axis_length
    return largest_region.minor_axis_length

    # The short axis is the minor axis length
    # short_axis_length = minor_axis_length


def NumberOfVoxels(image: sitk.Image):
    """
    Calculate the number of voxels in a binary mask image.

    Parameters:
    - image: SimpleITK image object representing the binary mask.
    """
    # Convert image to numpy array
    image_array = sitk.GetArrayFromImage(image)

    # Count the number of voxels in the binary mask
    num_voxels = (image_array > 0).sum()

    return num_voxels


def ActualVolume(image: sitk.Image):
    """
    Calculate the volume of a binary mask image.

    Parameters:
    - image_path: Path to the input binary mask image (e.g., NIfTI file).

    Returns:
    - volume: Volume of the binary mask in cubic millimeters.
    """
    # Calculate number of voxels
    num_voxels = NumberOfVoxels(image)

    # Calculate single voxel volume
    spacing = image.GetSpacing()  # (x, y, z) spacing in mm
    voxel_volume = spacing[0] * spacing[1] * spacing[2]  # in mm^3

    # Calculate total volume
    volume = num_voxels * voxel_volume  # in mm^3

    return volume


def get_measurements(image_path: Path):
    image = sitk.ReadImage(str(image_path))

    return {
        "Short_Axis_Length": ShortAxisDiameter(image),
        "Number_of_Voxels": NumberOfVoxels(image),
        "Actual_Volume (mm^3)": ActualVolume(image)
    }


@click.command()
@click.argument("image_dir", type=click.Path(exists=True, file_okay=False))
@click.option("--mask_type", type=click.Choice(['ground_truth', 'predicted']), default='predicted')
@click.option("--parallel", is_flag=True, help="Enable parallel processing")
@click.option("--n_jobs", type=int, default=-1, help="Number of parallel jobs to run (default: all available cores)")
def calc_image_props(image_dir: Path, mask_type: str, parallel: bool, n_jobs: int = -1):
    """
    Calculate the short axis length for all binary mask images in a directory
    and save the results to a csv file in the output directory.

    Parameters:
    - image_dir: Directory containing binary mask images.
    - output_dir: Directory where the output text files will be saved.
    """
    output_dir = Path(image_dir).parent

    mask_measurements = {}

    if parallel:
        image_paths = sorted(Path(image_dir).rglob("*GTVp*.nii*"))
        results = Parallel(n_jobs=n_jobs)(delayed(get_measurements)(image_path) for image_path in tqdm.tqdm(image_paths, desc="Processing masks in parallel"))
        mask_measurements = {image_path: result for image_path, result in zip(image_paths, results)}
    else:
        for image_path in tqdm.tqdm(sorted(Path(image_dir).rglob("*GTVp*.nii*")), desc="Processing masks"):
            mask_measurements[image_path] = get_measurements(image_path)

    # Save all mask measurements to a CSV file
    mask_measurements_df = pd.DataFrame.from_dict(mask_measurements, orient='index', columns=[
        "Short_Axis_Length", "Number_of_Voxels", "Actual_Volume (mm^3)"])
    mask_measurements_df.to_csv(output_dir / f"mask_measurements_{mask_type}.csv")

    return mask_measurements_df


if __name__ == "__main__":
    calc_image_props()  # pylint: disable=no-value-for-parameter
