import pandas as pd
import pydicom
import click
from damply import dirs
from pathlib import Path
from readii.io.loaders import loadImageDatasetConfig

def find_graphic_data(dataset: pydicom.Dataset,
                      file_path: str,
                      current_references: dict | None = None, 
                      graphic_data_list: list | None = None, 
                      ) -> list[dict]:
    """
    Recursively searches through a DICOM dataset to find and extract graphic data and related references.
    """
    if current_references is None:
        current_references = {}

    if graphic_data_list is None:
        graphic_data_list = []

    for elem in dataset:
        if elem.VR == 'SQ':  # Check for sequences
            for item in elem.value:
                find_graphic_data(item, file_path, current_references, graphic_data_list)
        elif elem.tag == (0x0070, 0x0022):  # Graphic Data tag
            graphic_type = current_references.get('GraphicType', 'Unknown')
            content_description = current_references.get('ContentDescription', 'Unknown')
            if graphic_type == 'Unknown' and content_description != 'Unknown':
                graphic_type = content_description
            graphic_info = {
                'SR Filename': file_path,
                'Graphic Type': graphic_type,
                'Reference Tag': current_references.get('ReferenceTag', 'Unknown'),
                'Reference': current_references.get('Reference', 'Unknown'),
                'Series Instance UID': current_references.get('SeriesInstanceUID', 'Unknown'),
                'x1': elem.value[0] if len(elem.value) > 0 else None,
                'y1': elem.value[1] if len(elem.value) > 1 else None,
                'x2': elem.value[2] if len(elem.value) > 2 else None,
                'y2': elem.value[3] if len(elem.value) > 3 else None
            }
            graphic_data_list.append(graphic_info)
        elif elem.tag == (0x0008, 0x1140):  # Referenced Image Sequence
            for ref_item in elem.value:
                sop_instance_uid = ref_item.ReferencedSOPInstanceUID
                current_references['ReferenceTag'] = 'ReferencedSOPInstanceUID'
                current_references['Reference'] = sop_instance_uid
                if 'ReferencedFrameNumber' in ref_item:
                    current_references['ReferencedFrameNumber'] = ref_item.ReferencedFrameNumber
                if 'SeriesInstanceUID' in ref_item:
                    current_references['SeriesInstanceUID'] = ref_item.SeriesInstanceUID
        elif elem.tag == (0x3006, 0x0024):  # Referenced Frame of Reference UID
            current_references['ReferenceTag'] = 'FrameOfReferenceUID'
            current_references['Reference'] = elem.value
        elif elem.tag == (0x0008, 0x1155):  # Referenced SOP Instance UID
            current_references['ReferenceTag'] = 'ReferencedSOPInstanceUID'
            current_references['Reference'] = elem.value
        elif elem.tag == (0x0008, 0x1160):  # Referenced Frame Number
            current_references['ReferenceTag'] = 'ReferencedFrameNumber'
            current_references['Reference'] = elem.value
        elif elem.tag == (0x0020, 0x000E):  # Series Instance UID
            current_references['SeriesInstanceUID'] = elem.value
        elif elem.tag == (0x0070, 0x0023):  # Graphic Type
            current_references['GraphicType'] = elem.value
        elif elem.tag == (0x0070, 0x0084):  # Content Description
            current_references['ContentDescription'] = elem.value


def extract_graphic_data_with_references(file_path: str) -> list[dict]:
    """
    Extract all graphic data and their references from a DICOM SR file.
    """
    ds = pydicom.dcmread(file_path)
    graphic_data_list = []
    find_graphic_data(ds, file_path, graphic_data_list = graphic_data_list)
    return graphic_data_list


@click.command
@click.option('--dataset', help='Name of dataset to process. Must have config file in config directory.')
@click.option('--overwrite', help='Whether to overwrite existing readii image files', default=False)
def scrape_annotation_dataset(dataset:str,
                              overwrite:bool = False) -> pd.DataFrame:
    """
    Scrape all DICOM SR files for a given dataset to make an index file.
    """
    # Load in configuration settings for dataset
    config = loadImageDatasetConfig(dataset, dirs.CONFIG)

    # Set up output path
    out_path = dirs.PROCDATA / f"{config['DATA_SOURCE']}_{config['DATASET_NAME']}" / "annotations" / f"{config['DATASET_NAME']}_graphic_data.csv"

    # Check if output already exists
    if out_path.exists() and not overwrite:
        print(f"{dataset} annotations have already had their graphic data scraped. See {out_path}. Set overwrite to True to force scrape. Loading existing file.")
        dataset_graphic_df = pd.read_csv(out_path)

    else:
        # Create output path directories
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Set up path to annotation files
        annotation_dir_path = dirs.RAWDATA / f"{config['DATA_SOURCE']}_{config['DATASET_NAME']}" / "images" / "annotations" / f"{config['DATASET_NAME']}"
        # Initialize empty list to store graphic data in for each file
        dataset_graphic_data = []
        # Iterate over the annotation files, scrape their graphic data, and save the list in the global list. Should end up as a list of dictionaries.
        for file_path in annotation_dir_path.glob("*.dcm"):
            dataset_graphic_data += (extract_graphic_data_with_references(file_path))
        
        # Convert list of dictionaries to a dataframe
        dataset_graphic_df = pd.DataFrame(dataset_graphic_data)
        # Save out dataframe
        dataset_graphic_df.to_csv(out_path, index=False)

    return dataset_graphic_df
    
        

if __name__ == '__main__':
    # Example usage
    # dicom_file_path = "data/rawdata/TCIA_CPTAC-HNSCC/images/annotations/CPTAC-HNSCC/animated.peafowl-CPTAC-HNSCC-C3N-00828-1.3.6.1.4.1.14519.5.2.1.3320.3273.324519151134894128406738207272.dcm"
    # graphic_data_df = extract_graphic_data_with_references(dicom_file_path)
    # print(graphic_data_df)

    scrape_annotation_dataset()
