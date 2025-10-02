#!/bin/bash
# Using a configuration file (.yaml) and imgtools, process images into NIFTIs

### Keyboard inputs ###
# Help information 
help_funct()
{
    echo ""
    echo "Usage of $0: -a DISEASE_SITE -b DATASET_NAME -c YAML_PATH -d EXIST_FILE"
    echo -e "\t-a DISEASE_SITE: Area where the disease is, which is the folder name above the dataset name folder"
    echo -e "\t-b DATASET_NAME: Full name of dataset (e.g. TCIA_CPTAC-CCRCC)"
    echo -e "\t-c YAML_PATH: Path to .yaml configuration file" 
    echo -e "\t-d EXIST_FILE: med-imagetools parameter determining what to do if files exist already. Options are overwrite, skip, fail."
    exit 1
}
# Get necessary args 

while getopts "a:b:c:d:?" opt
do 
    case "$opt" in
        a ) DISEASE_SITE="$OPTARG" ;;
        b ) DATASET_NAME="$OPTARG" ;;
        c ) YAML_PATH="$OPTARG" ;;
        d ) EXIST_FILE="$OPTARG" ;;
        ? ) help_funct ;; #If inputted parameter does not exist, show information in help function
    esac
done 

# Check if parameters are empty 
if [ -z "$DISEASE_SITE" ] || [ -z "$DATASET_NAME" ] || [ -z "$YAML_PATH" ] || [ -z "$EXIST_FILE" ] ;
then 
    echo "Some necessary parameters were not defined/left empty Please see below for help";
    help_funct
fi 

### Executing the image processing ###
# Read info from configuration file
# -r needed else the strings will include double quotes around them
IMG_SPACING=$(yq -r '.img_settings.spacing' $YAML_PATH)
ROI_MATCH=$(yq -r '.file_settings.roi_match_substring' $YAML_PATH)
FILE_FORMAT=$(yq -r '.file_settings.file_format' $YAML_PATH)
MODALITIES=$(yq -r '.file_settings.modalities' $YAML_PATH)
UNIQUE_OUT_ID=$(yq -r '.file_settings.unique_out_id' $YAML_PATH)

# Separate dataset name into where it's from and what it's called (used to define output file name)
# split_dataname[0] is where it's from (e.g. PMCC, TCIA, etc.)
# split_dataname[1] is what the dataset is called (e.g. CPTAC-CCRCC, NSCLC-Radiogenomics, etc.)
IFS='_' read -ra split_dataname <<< $DATASET_NAME
DATANAME=${split_dataname[1]}

# Execute imgtools to process dcm to nifti with given parameters
pixi run imgtools autopipeline \
$RAWDATA/$DISEASE_SITE/$DATASET_NAME/images \
$PROCDATA/$DISEASE_SITE/$DATASET_NAME/images/mit\_$DATANAME\_$UNIQUE_OUT_ID \
--modalities $MODALITIES \
--filename-format $FILE_FORMAT \
--existing-file-mode $EXIST_FILE \
--roi-match-map $ROI_MATCH \
--spacing $IMG_SPACING \
--update-crawl