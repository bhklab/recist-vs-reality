# Usage Guide

## Project Configuration

Each dataset needs a configuration file with the following settings filled in

```
DATA_SOURCE: ""    # where the data came from, will be used for data organization
DATASET_NAME: ""   # the name of the dataset , will be use for data organization

### MED-IMAGETOOLS settings
MIT:
    MODALITIES:                 # Modalities to process with autopipeline
        image: CT
        mask: RTSTRUCT     
    ROI_STRATEGY: MERGE         # How to handle multiple ROI matches 
    ROI_MATCH_MAP:              # Matching map for ROIs in dataset (use if you only want to process some of the masks in a segmentation)
        KEY:ROI_NAME            # NOTE: there can be no spaces in KEY:ROI_NAME
```

## Data Setup
The following sections describe how to set up the data you wish to process with this pipeline following the BHKLab Data Management Protocol (DMP). This will ensure data remains separate from the project directory and accessible to other users.

### Raw Data
Set up a separate main data directory outside of the project directory. We'll call this `Datasets`.

In `Datasets`, set up a directory for the dataset you wish to process as follows:

```bash
Datasets
|---- {DATASET_SOURCE}_{DATASET_NAME}
      |-- clinical
      |   `-- {Clinical Data File}.csv OR {Clinical Data File}.xlsx
      `-- images
          |-- {DATASET_NAME}
          |   |-- {PatientID}
          |   |   `-- {StudyUID}
          |   |       |-- {Image DICOM directory}
          |   |       |   |-- 1-01.dcm
          |   |       |   |-- ...
          |   |       |   |-- 1-N.dcm
          |   |       |-- {Mask DICOM directory}
          |   |       |   `-- 1-01.dcm
          |   |-- {PatientID}
          |   |-- ...
          |   `-- {PatientID}
          `-- annotations
              `-- {DATASET_NAME}
                  |-- DICOM-SR_annotation_file.dcm
                  |-- DICOM-SR_annotation_file.dcm
                  `-- DICOM-SR_annotation_file.dcm
```
Image directory structure may vary depending on the source. This example is based on the structure setup by TCIA when downloading with a manifest file.

!!! note "BHKLab DMP Setup"
    If using the BHKLab DMP, the `Datasets` directory will be structured with `rawdata/{DiseaseRegion}/{DATASET_SOURCE}_{DATASET_NAME}`. In the next step, you can create the symbolic link starting from `{DATASET_SOURCE}_{DATASET_NAME}`.

Once this data directory is setup, run the following in a terminal from the main directory of the project.

```bash
ln -s /path/to/Datasets/{DATASET_SOURCE}_{DATASET_NAME} data/rawdata
```

This will create a symbolic link to your dataset in the `Datasets` directory. 

You can confirm this worked by running:
```bash
ls -l data/rawdata
```
and you should see,
```bash
total 5
-rw-rw-r-- 1 bhkuser root 1395 Jun  4 15:21 README.md
lrwxrwxrwx 1 bhkuser root   80 Jun  4 15:46 {DATASET_SOURCE}_{DATASET_NAME} -> /path/to/Datasets/{DATASET_SOURCE}_{DATASET_NAME}
```

Now, document the dataset you've added on the [Data Sources](data_sources.md) page following the provided template.

---

### Processed Data
If you wish to use the BHKLab DMP strategy, follow the process below.

Create a processed data directory for your dataset in the external `Datasets` directory as follows:
```bash
mkdir /path/to/Datasets/procdata/{DiseaseRegion}/{DATASET_SOURCE}_{DATASET_NAME}
```
**Note**: the Disease Region must match what is in the `rawdata` path to the dataset.

Now you can create a symbolic link to this directory in the project directory, like we did for the raw data.
```bash
ln -s /path/to/Datasets/procdata/{DiseaseRegion}/{DATASET_SOURCE}_{DATASET_NAME} data/procdata
```

### Results
TODO:: describe results directory setup


## Running Your Analysis

### 1. Running Med-ImageTools

TODO:: discuss using pixi tasks to run the analysis
