# Developer Notes - Kaitlyn 

## Annotation-Imaging-Segmentation Matching 

### Additional Function Documentation (Last Updated [2025-11-13])

#### get_ann_measurements() 

The DICOM data extracted through `pydicom.dcmread` is in a similar structure to nested dictionaries. Basic information (PatientID, annotation ID, etc.) can be found from the outer dictionary. The reference imaging information is found within the Current Requested Procedure Evidence Sequence one level down from the outer dictionary. 

The annotation information is found within nested content sequences, with the outermost content sequence (I've called the parent content sequence) usually containing 5 items with information about language, country of origin, etc., the file itself, the procedure(s) used, and measurements. Measurements are usually the last item in the parent sequence (at index 4). Each measurement is kept in another content sequence. The number of these content sequences are determined by the number of measurements taken. Must loop through all of them to get all annotation measurements. The content sequence each measurement is found in usually has a length of 5. The long axis measurement information is found in the 3rd item (index 2) and the short axis measurement information is found in the 4th item (index 3). 

Once inside either the long axis or short axis content sequence, you can get the following corresponding information: 
* Measurement Type: Found by searching the Concept Name Code Sequence for the first item and then obtaining the code meaning. The value of the measurement type should be "Long Axis" or "Short Axis" depending on what is being looked for. 
* Measurement Unit: Founed by searching the Measured Value Sequence for the first item, then searching for the Measurement Units Code Sequence first item, and then accessing the Code Value. This is usually "mm". 
* Annotation Measurement: Found by searching the Measured Value Sequence for the first item (same as above) and accessing the Numeric Value 
* Referenced SOPUID: Denotes the slice that the measurement was taken on (paired long and short axis measurements should be on the same slice). Found by entering the first item of two additional content sequences and the Referenced SOP Sequence and then accessing the Referenced SOP Instance UID value. This number should correspond to an instance value in the matching image (listed for each image in the crawl_db.json file). 
* Axis Coordinates: Stored as x1/y1/x2/y2 in the DICOM data, these are the points that define the annotation drawn. Found by searching the first item of one additional content sequence and looping through the values of Graphic Data, which should always have a length of 4. 

#### get_rtstruct_SOPUIDs(): 

This function is currently necessary due to an issue with med-imagetools not extracting the slice IDs associated with RTSTRUCT segmentation files specifially. This type of file will not have any ReferencedSOPUIDs listed in the crawl_db.json function, unlike the SEG files. These SOPUIDs are necessary to match the annotation with the segmentation. 

Similarly to `get_ann_measurements()`, the referenced slice IDs are found in a nested dictionary of sequences, each with only one item in them. The Referenced Frame of Reference Sequence contains information related to the slices in the segmentation file, the frame of reference, and the image the segmentation pertains to. **The list of slices may only contain slices with actual segmentation or it could contain all slices within the referenced image (including where the binary mask is all 0 for a slice)**. 

Once inside the Referenced Frame of Reference Sequence, the RT Referenced Study Sequence of length 1 will need to be entered. This includes all information in the previous sequence except for the frame of reference information. At this level, you can access the study instance UID and the RT Referenced Series Sequence with length 1. This contains the information related to the imaging series instance UID and the segmentation contour. The contour information is stored in a Contour Image Sequence whose length is determined by however many slices the segmentation is. This sequence must be looped through to find each slice ID (found with the ReferencedSOPInstanceUID tag). 

#### get_slice_num(): 

The instance name of the slice being referenced usually is in the form of "subseries_num-slice_num.dcm" (e.g. 1-001.dcm). At this stage, only the slice that was referenced in an annotation will be passed into this function. To get the slice number, the instance name first gets split by "." to remove the extention and then gets split by "-" to remove the subseries. 

The nifti files are ordered backwards compared to the DICOM slices and they differ by having a starting index of 0 instead of 1. To make sure the correct nifti slice is selected for the measurement and segmentation, the found slice number gets subtracted from the total number of slices to get the position of the desired nifti slice. 

## Preprocessing for MedSAM2-RECIST 
### Additional Function Documentation (Last Updated [2025-11-13])

#### create_npzs()
The MedSam2-RECIST model requires that at there at least be 5 points to be sampled from the RECIST annotation. That means that if there are less than 5 pixels involved in the annotation, an error will be produced and the run will stop. Therefore, a theshold of 5 pixels is used to see if the data should be exported to .npz format. 