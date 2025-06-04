# Data Sources

## External Data Sources

### CPTAC-HNSCC
- **Name**: The Clinical Proteomic Tumor Analysis Consortium Head and Neck Squamous Cell Carcinoma Collection
- **Version/Date**: Version 14 or 15?
- **URL**: <https://www.cancerimagingarchive.net/collection/cptac-hnscc/>
- **Access Method**: NBIA Data Retriever
- **Access Date**: 2024-03-29
- **Data Format**: DICOM-CT, DICOM-RTSTRUCT
- **Citation**: National Cancer Institute Clinical Proteomic Tumor Analysis Consortium (CPTAC). (2018). The Clinical Proteomic Tumor Analysis Consortium Head and Neck Squamous Cell Carcinoma Collection (CPTAC-HNSCC) (Version 19) [dataset]. The Cancer Imaging Archive. https://doi.org/10.7937/k9/tcia.2018.uw45nh81
- **License**: [TCIA Restricted License](https://wiki.cancerimagingarchive.net/download/attachments/4556915/TCIA%20Restricted%20License%2020220519.pdf?api=v2)
- **Data Types**: 
    - Images: CT, RTSTRUCT
- **Sample Size Selected**: TBD
- **PatientIDs**: 
    * C3N-00297
    * C3N-00498
    * C3N-00828
    * C3N-01620
    * C3N-01752
    * C3N-01754
    * C3N-01757
    * C3N-01943
    * C3N-01948

### Crowds-Cure-2018
- **Name**: Crowds Cure Cancer: Data collected at the RSNA 2018 annual meeting
- **Version/Date**: Version 1: Updated 2019/05/30
- **URL**: <https://www.cancerimagingarchive.net/analysis-result/crowds-cure-2018/>
- **Access Method**: NBIA Data Retriever
- **Access Date**: 2025-03-26
- **Data Format**: DICOM-SR
- **Citation**: Urban, T., Ziegler, E., Pieper, S., Kirby, J., Rukas, D., Beardmore, B., Somarouthu, B., Ozkan, E., Lelis, G., Fevrier-Sullivan, B., Nandekar, S., Beers, A., Jaffe, C., Freymann, J., Clunie, D., Harris, G. J., & Kalpathy-Cramer, J. (2019). Crowds Cure Cancer: Crowdsourced data collected at the RSNA 2018 annual meeting [Data set]. The Cancer Imaging Archive. https://doi.org/10.7937/TCIA.2019.yk0gm1eb
- **License**: [CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)
- **Data Types**: 
    - Annotations: Structured Reports containing 2D RECIST measurements
- **Sample Size Selected**: TBD


## Overview

This section should document all data sources used in your project.
Proper documentation ensures reproducibility and helps others
understand your research methodology.

## How to Document Your Data

For each data source, include the following information:

### 1. External Data Sources

- **Name**: Official name of the dataset
- **Version/Date**: Version number or access date
- **URL**: Link to the data source
- **Access Method**: How the data was obtained (direct download, API, etc.)
- **Access Date**: When the data was accessed/retrieved
- **Data Format**: Format of the data (FASTQ, DICOM, CSV, etc.)
- **Citation**: Proper academic citation if applicable
- **License**: Usage restrictions and attribution requirements

Example:

```markdown
## TCGA RNA-Seq Data

- **Name**: The Cancer Genome Atlas RNA-Seq Data
- **Version**: Data release 28.0 - March 2021
- **URL**: https://portal.gdc.cancer.gov/
- **Access Method**: GDC Data Transfer Tool
- **Access Date**: 2021-03-15
- **Citation**: The Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. Nature, 490(7418), 61-70.
- **License**: [NIH Genomic Data Sharing Policy](https://sharing.nih.gov/genomic-data-sharing-policy)
```

### 2. Internal/Generated Data

- **Name**: Descriptive name of the dataset
- **Creation Date**: When the data was generated
- **Creation Method**: Brief description of how the data was created
- **Input Data**: What source data was used
- **Processing Scripts**: References to scripts/Github Repo used to generate this data

Example:

```markdown
## Processed RNA-Seq Data
- **Name**: Processed RNA-Seq Data for TCGA-BRCA
- **Creation Date**: 2021-04-01
- **Creation Method**: Processed using kallisto and DESeq2
- **Input Data**: FASTQ Data obtained from the SRA database
- **Processing Scripts**: [GitHub Repo](https://github.com/tcga-brca-rnaseq)
```

### 3. Data Dictionary

For complex datasets, include a data dictionary that explains:

| Column Name | Data Type | Description | Units | Possible Values |
|-------------|-----------|-------------|-------|-----------------|
| patient_id  | string    | Unique patient identifier | N/A | TCGA-XX-XXXX format |
| age         | integer   | Patient age at diagnosis | years | 18-100 |
| expression  | float     | Gene expression value | TPM | Any positive value |

## Best Practices

- Store raw data in `data/rawdata/` and never modify it
- Store processed data in `data/procdata/` and all code used to generate it should be in `workflow/scripts/`
- Document all processing steps
- Track data provenance (where data came from and how it was modified)
- Respect data usage agreements and licenses!
    This is especially important for data that should not be shared publicly
