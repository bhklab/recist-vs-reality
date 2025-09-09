import pandas as pd 
import numpy as np 
import click

from damply import dirs
from pathlib import Path

def get_tumour_measurements(feature_df: pd.DataFrame, 
                            diameters: str): 
    '''
    Retrieves all patient and file information, diagnostic information related to imaging spacing and size, and volume and diameter measurements. 

    Parameters 
    ----------
    feature_df: pd.DataFrame
        The pyRadiomics feature dataframe you would like to obtain a subset of data from 
    diameters: str 
        Decides which set of diameter measurements is used for tumour volume analysis (see function calc_gordon_tva)
            max_diam --> Maximum 2D diameters 
            axis_len --> Axis lengths 
    
    Returns
    ----------
    pat_diag_feat_df: pd.DataFrame 
        A dataframe containing only the specified subset of information and features
    '''

    def get_pat_and_file_info(feat_df): 
        '''
        Obtains all patient and file information ahead of the diagnostics information in the feature file. 

        Parameters
        ----------
        feat_df: pd.DataFrame
            The pyRadiomics feature dataframe you would like to obtain a subset of data from
        
        Returns
        ----------
        pat_file_info: pd.DataFrame 
            A subset of feature_df containing only the patient and file information along with their associated indices
        '''
        diag_info_start = feat_df.columns.get_loc("diagnostics_Versions_PyRadiomics")
        pat_file_info = feat_df.iloc[:, 0:diag_info_start]

        return pat_file_info

    def calc_gordon_tva(feat_df, diams): 
        '''
        Calculates the tumour volume analysis (TVA) outlined in Gordon et al: https://doi.org/10.1016/j.ejso.2025.109578
        TVA should be calculated as 4/3*pi*(MajorAxisLength)/2*(MinorAxisLength)/2*(LeastAxisLength)/2
        
        Parameters
        ----------
        feat_df: pd.DataFrame
            The pyRadiomics feature dataframe you would like to obtain a subset of data from
        diams: str 
            Decides which set of diameter measurements is used
                max_diam --> Maximum 2D diameters (OLD, DO NOT USE)
                axis_len --> Axis lengths 
        
        Returns 
        ----------
        feat_w_tva_df: pd.DataFrame 
            The original dataframe with the addition of a column containing the TVA values for all patients
        '''
        if diams == "max_diam":
            #This calculation is not for an encompassing ellipsoid. Will likely get removed later
            feat_df["original_shape_gordon_tva"] = 4/3*np.pi*feat_df["original_shape_Maximum2DDiameterRow"]*feat_df["original_shape_Maximum2DDiameterColumn"]*feat_df["original_shape_Maximum2DDiameterSlice"]
        elif diams == "axis_len":
            feat_df["original_shape_gordon_tva"] = 4/3*np.pi*feat_df["original_shape_MajorAxisLength"]/2*feat_df["original_shape_MinorAxisLength"]/2*feat_df["original_shape_LeastAxisLength"]/2 #Divide to get semi-axis lengths
        else: 
            raise ValueError("Diameter setting not accepted. Please choose either max_diam or axis_len.")
        feat_w_tva_df = feat_df.copy(deep = True) 

        return feat_w_tva_df

    feat_df_w_tva = calc_gordon_tva(feature_df, diams = diameters)
    patient_file_info = get_pat_and_file_info(feature_df)

    diagnostic_labels = ["diagnostics_Image-original_Spacing", 
                         "diagnostics_Image-interpolated_Spacing",
                         "diagnostics_Image-original_Size", 
                         "diagnostics_Image-interpolated_Size"
                        ]
    
    feature_labels = ["original_shape_MeshVolume", 
                      "original_shape_VoxelVolume", 
                      "original_shape_Maximum3DDiameter", 
                      "original_shape_Maximum2DDiameterColumn", 
                      "original_shape_Maximum2DDiameterRow", 
                      "original_shape_Maximum2DDiameterSlice", 
                      "original_shape_MajorAxisLength",
                      "original_shape_MinorAxisLength",
                      "original_shape_LeastAxisLength",
                      "original_shape_gordon_tva"
                     ]
    
    labels_of_interest = diagnostic_labels + feature_labels

    diag_feat_df = feat_df_w_tva[labels_of_interest]

    pat_diag_feat_df = pd.concat([patient_file_info, diag_feat_df], axis = 1)

    return pat_diag_feat_df 

@click.command()
@click.option('--feat_path', help = 'Path to pyradiomics feature .csv file')
@click.option('--export_path', help = 'Path to export directory')
@click.option('--export_filename', help = 'Name of consolidated feature .csv file')
def create_rvr_measurements(feat_path: str, 
                            export_path: str,
                            export_filename: str): 
    '''
    Consolidates measurement features from pyradiomics feature file and exports it to a provided destination. 
    
    Parameters
    ----------
    feat_path: str
        Path to the pyradiomics feature file csv
    export_path: str
        Output file path
    export_filename: str
        Filename to be concatenated with the export path
    '''
    feature_data = pd.read_csv(feat_path)

    consol_feats = get_tumour_measurements(feature_data, diameters = "axis_len")

    export_path = Path(export_path)
    if not export_path.exists():
        Path(export_path).mkdir(parents=True, exist_ok = True)

    consol_feats.to_csv(export_path / export_filename, index = False)

if __name__ == '__main__':
    create_rvr_measurements()