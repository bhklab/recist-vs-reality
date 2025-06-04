import pydicom
import pandas as pd

def find_graphic_data(dataset, current_references, graphic_data_list, file_path):
    """
    Recursively searches through a DICOM dataset to find and extract graphic data and related references.
    """
    for elem in dataset:
        if elem.VR == 'SQ':  # Check for sequences
            for item in elem.value:
                find_graphic_data(item, current_references, graphic_data_list, file_path)
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

def extract_graphic_data_with_references(file_path):
    """
    Extract all graphic data and their references from a DICOM SR file.
    """
    ds = pydicom.dcmread(file_path)
    graphic_data_list = []
    find_graphic_data(ds, {}, graphic_data_list, file_path)
    return pd.DataFrame(graphic_data_list)

# Example usage
dicom_file_path = "/home/bhkuser3/RECIST-play/rawdata/structreports/OCT01-1123-COB4702_1.dcm"
graphic_data_df = extract_graphic_data_with_references(dicom_file_path)
