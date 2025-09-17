
# Developer Notes - Kaitlyn 

## Annotation-Imaging-Segmentation Matching 

### Additional Function Documentation (Last Updated [2025-09-17])

#### get_ann_measurements()

The DICOM data extracted through `pydicom.dcmread` is in a similar structure to nested dictionaries. Basic information (PatientID, annotation ID, etc.) can be found from the outer dictionary. The reference imaging information is found within the Current Requested Procedure Evidence Sequence one level down from the outer dictionary. 

The annotation information is found within nested content sequences, with the outermost content sequence (I've called the parent content sequence) usually containing 5 items with information about language, country of origin, etc., the file itself, the procedure(s) used, and measurements. Measurements are usually the last item in the parent sequence (at index 4). Each measurement is kept in another content sequence. The number of these content sequences are determined by the number of measurements taken. Must loop through all of them to get all annotation measurements. Each measurement is found in another content sequence, usually with a length of 5. The long axis measurement information is found in the 3rd item (index 2) and the short axis measurement information is found in the 4th item (index 3). 

Once inside either the long axis or short axis content sequence, you can get the following corresponding information: 
* Measurement Type: Found by searching the Concept Name Code Sequence for the first item and then obtaining the code meaning. The value of the measurement type should be "Long Axis" or "Short Axis" depending on what is being looked for. 
* Measurement Unit: Founed by searching the Measured Value Sequence for the first item, then searching for the Measurement Units Code Sequence first item, and then accessing the Code Value. This is usually "mm". 
* Annotation Measurement: Found by searching the Measured Value Sequence for the first item (same as above) and accessing the Numeric Value 
* Referenced SOPUID: Denotes the slice that the measurement was taken on (paired long and short axis measurements should be on the same slice). Found by entering the first item of two additional content sequences and the Referenced SOP Sequence and then accessing the Referenced SOP Instance UID value. This number should correspond to an instance value in the matching image (listed for each image in the crawl_db.json file). 
* Axis Coordinates: Stored as x1/y1/x2/y2 in the DICOM data, these are the points that define the annotation drawn. Found by searching the first item of one additional content sequence and looping through the values of Graphic Data, which should always have a length of 4. 