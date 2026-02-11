from pathlib import Path

def get_source_dataset(source_id:str) -> str:
    # Get dataset name based on source_id value (e.g.KiTS23_case_00001)
    if "KiTS" in source_id: 
        dataset = "KiTS"
    elif "LIDC-IDRI" in source_id: 
        dataset = "LIDC-IDRI"
    elif "LNDb" in source_id: 
        dataset = "LNDb" 
    elif "MSD_colon" in source_id: 
        dataset = "MSD Colon"
    elif "MSD_hepaticvessel" in source_id: 
        dataset = "MSD Hepatic\nVessel" 
    elif "MSD_liver" in source_id: 
        dataset = "MSD Liver" 
    elif "MSD_lung" in source_id: 
        dataset = "MSD Lung"
    elif "MSD_pancreas" in source_id: 
        dataset = "MSD Pancreas"
    elif "NSCLC-Radiogenomics" in source_id: 
        dataset = "NSCLC-\nRadiogenomics" 
    else: 
        dataset = 'unknown'
    return dataset


def add_source_labels(main_df, mapping_df):
    # Merge this name mapping into the main dataframe for med
    main_lbl_df = main_df.merge(mapping_df, left_on="sample_id", right_on="File Name", how="left")
    # Drop duplicate column File Name, is same as sample_id
    main_lbl_df = main_lbl_df.drop(labels="File Name", axis=1)
    # Rename Source column to source_id for continuity with sample_id
    main_lbl_df = main_lbl_df.rename(columns={'Source': "source_id"})

    # Add source_dataset column to main dataframe based on source_id
    main_lbl_df['source_dataset'] = main_lbl_df["source_id"].apply(get_source_dataset)

    return main_lbl_df


def add_sample_id_column(main_df, 
                         image_path_col="image_path"):
    # Get sample_id from image_path column, which is the parent directory name of the image file
    def get_sample_id(path):
        return Path(path).parent.stem
    
    main_df['sample_id'] = main_df[image_path_col].apply(get_sample_id)
    return main_df