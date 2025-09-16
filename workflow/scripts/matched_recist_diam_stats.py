import pandas as pd 
import logging
import matplotlib.pyplot as plt

from scipy import stats
from pathlib import Path 
from damply import dirs

def get_nifti_locs(nifti_idx_path: Path, 
                   img_seg_info: pd.DataFrame): 
    '''
    Gets the file locations of the nifti files for CT image and associated segmentation. 
    
    Parameters
    ----------
    nifti_idx_path: Path
        Path to index file created by med-imagetools for the nifti versions of the CT and segmentation files  
    img_seg_info: pd.DataFrame 
        A one-row dataframe containing the current image and segmentation information

    Returns 
    ----------
    img_nifti_path: Path 
        Path to the CT imaging nifti file 
    seg_nifti_path: Path 
        Path to the associated segmentation nifti file 
    '''
    no_nifti = False
    nifti_path = Path("/".join(str(nifti_idx_path).split("/")[:-1]))
    nifti_idx_df = pd.read_csv(nifti_idx_path)

    img_instUID = img_seg_info["ImgSeriesInstanceUID"].values[0]
    seg_instUID = img_seg_info["SegSeriesInstanceUID"].values[0]

    img_nifti = nifti_idx_df[nifti_idx_df["SeriesInstanceUID"] == img_instUID]["filepath"]
    seg_nifti = nifti_idx_df[nifti_idx_df["SeriesInstanceUID"] == seg_instUID]["filepath"]

    if img_nifti.empty: 
        logger.info("Imaging SeriesInstanceUID: %s does not have a matching nifti file.")
        no_nifti = True
    elif seg_nifti.empty: 
        logger.info("Segmentation SeriesInstanceUID: %s does not have a matching nifti file.")
        no_nifti = True

    if no_nifti: 
        return "", ""
    img_nifti_path = nifti_path / img_nifti.values[0]
    seg_niftis = nifti_path / seg_nifti

    for row in seg_niftis.items(): 
        path = str(row[1])
        if "GTV" in path: 
            seg_nifti_path = Path(path)
    
    return img_nifti_path, seg_nifti_path 

def match_pyrad_to_ann_seg(nifti_file_path: Path, 
                           pyrad_data: pd.DataFrame, 
                           ann_seg_match_data: pd.DataFrame): 
    '''
    Match the measurements taken with PyRadiomics with the appropriate annotation and segmentation match. 

    Parameters
    ----------
    nifti_file_path: Path, 
        Path to the nifti index.csv file created by med-imagetools 
    pyrad_data: pd.DataFrame, 
        Contains the consolidated measurement data from PyRadiomics 
    ann_seg_match_data: pd.DataFrame, 
        Contains the matched annotation, image, and segmentation data 
    
    Returns
    ----------
    ann_seg_pyrad_data: pd.DataFrame, 
        Contains the matched annotation, image, segmentation, and PyRadiomics measurement data. Each annotation has 
        it's own row. 
    '''
    ann_seg_pyrad_data = pd.DataFrame()

    for idx, row in ann_seg_match_data.iterrows(): 
        curr_matched = row.to_frame().transpose().reset_index(drop = True) #Gets pd.Series data back into a dataframe with correct columns
    
        if len(curr_matched["SegSeriesInstanceUIDs"].values) == 1:
            #Get current segmentation ID since it's stored in an array
            curr_matched["SegSeriesInstanceUID"] = curr_matched["SegSeriesInstanceUIDs"].values[0]
        else: 
            #Will enter here if there were two segmentations of the same tumour that made it through the 
            #matching process. 
            #This shouldn't happen for the datasets processed so far (CCRCC, NSCLC-Radiogenomics, PDA)
            logger.debug("Current matched row has more than one segmentation for this tumour.") 
        
        _, seg_nifti_path = get_nifti_locs(nifti_file_path, curr_matched)

        #Get folder name of where the nifti files are located
        #This should be the start of the nifti segmentation path listed in the PyRadiomics feature .csv under the header "Mask"
        #Can be removed later if NSCLC-Radiogenomics nifti paths match
        nifti_folder = str(nifti_file_path).split("/")[-2]

        if "NSCLC-Radiogenomics" in nifti_folder: 
            #Additional processing of path since the features were extracted on a previous (but equivalent) nifti file
            #of a different name
            #Should be removed later if the features ever get re-extracted and the file paths reflect the current 
            #nifti files in the dataset folder
            seg_nifti_path = "/".join(str(seg_nifti_path).split("/")[-4:-1]) / "GTV.nii.gz"
        else: 
            #Get path that will match form of the PyRadiomics "Mask" path in the feature file
            seg_nifti_path = "/".join(str(seg_nifti_path).split("/")[-4:])
        
        #Check if the current segmentation was used in feature extraction. If so, match it with the appropriate annotation,
        #imaging, and segmentation data for further processing
        if pyrad_data["Mask"].isin([str(seg_nifti_path)]).any(): 
            curr_pyrad = pyrad_data[pyrad_data["Mask"] == str(seg_nifti_path)].reset_index(drop = True) 
            curr_match_pyrad = curr_matched.join(curr_pyrad)

            if ann_seg_pyrad_data.empty: 
                ann_seg_pyrad_data = curr_match_pyrad 
            else: 
                ann_seg_pyrad_data = pd.concat([ann_seg_pyrad_data, curr_match_pyrad]).reset_index(drop = True)

    return ann_seg_pyrad_data 

def calc_ttest_stats(ann_seg_pyrad_df: pd.DataFrame, 
                     feat_names: list, 
                     epsilons: list): 
    '''
    Calculates the p-values for a two-sample t-test and two one-sided t-tests (TOST) comparing the long axis annotation
    measurements and the desired PyRadiomics measurements.

    Parameters
    ----------
    ann_seg_pyrad_df: pd.DataFrame 
        Contains the matched annotation, image, segmentation, and PyRadiomics measurement data. Each annotation has 
        it's own row
    feat_names: list 
        The PyRadiomics features you would like to test against the long axis annotation measurements. Must be as they
        appear in the PyRadiomics feature file
    epsilons: list
        Specifies the equivalence margins that you would like to test in TOST

    Returns
    ----------
    ann_pyrad_stats: pd.DataFrame 
        Contains the statistical test results, with each row representing a PyRadiomics feature tested against the 
        long axis length
    '''

    ann_pyrad_stats = pd.DataFrame()

    for feature in feat_names: 
        curr_stats = pd.DataFrame()
        curr_stats["FeatureName"] = [feature]
        
        #Two sample t-test 
        two_samp_t = stats.ttest_ind(ann_seg_pyrad_df["AnnLongAxisLength"].astype(float), ann_seg_pyrad_df[feature])
        curr_stats["TwoSampleTTest_stat"] = [two_samp_t[0]]
        curr_stats["TwoSampleTTest_p"] = [two_samp_t[1]]
        # curr_stats["TwoSampleTTest_df"] = two_samp_t[2] #For some reason this gives an index out of range error? Degrees of freedom not a returnable value?

        #TOST
        for e in epsilons: 
            _, p_greater = stats.ttest_ind(ann_seg_pyrad_df["AnnLongAxisLength"].astype(float) + e, 
                                           ann_seg_pyrad_df[feature], 
                                           alternative = "greater")
            _, p_lesser = stats.ttest_ind(ann_seg_pyrad_df["AnnLongAxisLength"].astype(float) - e, 
                                          ann_seg_pyrad_df[feature], 
                                          alternative = "less")
            pval = max(p_greater, p_lesser) 

            header = "TOST_e" + str(e) + "_p" #Create header name for TOST p-value based on the epsilon value used
            curr_stats[header] = [pval]
        
        if ann_pyrad_stats.empty: 
            ann_pyrad_stats = curr_stats
        else: 
            ann_pyrad_stats = pd.concat([ann_pyrad_stats, curr_stats]).reset_index(drop=True)
    
    return ann_pyrad_stats

if __name__ == '__main__': 
    dataset = "TCIA_CPTAC-PDA" #Only used to help distinguish between log files currently

    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=dirs.LOGS / (dataset + "_recist_diam_stats.log"), encoding='utf-8', level=logging.DEBUG)

    nifti_idx_path = dirs.PROCDATA / "Abdomen/TCIA_CPTAC-PDA/images/mit_CPTAC-PDA/mit_CPTAC-PDA_index.csv"
    pyrad_measure_path = dirs.PROCDATA / "Abdomen/TCIA_CPTAC-PDA/features/rvr_measurements/pyradiomics_measurement_subset_axislen.csv"
    matched_data_path = dirs.PROCDATA / "Abdomen/TCIA_CPTAC-PDA/metadata/annotation_seg_matching/matching_ann_to_seg.csv"
    out_path = dirs.PROCDATA / "Abdomen/TCIA_CPTAC-PDA/metadata/recist_diameter_stats"

    pyrad_df = pd.read_csv(pyrad_measure_path) 
    matched_data_df = pd.read_csv(matched_data_path)

    match_pyrad_ann_df = match_pyrad_to_ann_seg(nifti_file_path = nifti_idx_path, 
                                                pyrad_data = pyrad_df, 
                                                ann_seg_match_data = matched_data_df)
    
    eps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #Trying different epsilons until we decide on an acceptable equivalence margin
    feature_names = ["original_shape_Maximum3DDiameter", 
                     "original_shape_Maximum2DDiameterColumn", 
                     "original_shape_Maximum2DDiameterRow", 
                     "original_shape_Maximum2DDiameterSlice", 
                     "original_shape_MajorAxisLength", 
                     "original_shape_MinorAxisLength", 
                     "original_shape_LeastAxisLength"]
    
    annotation_pyrad_stats = calc_ttest_stats(ann_seg_pyrad_df = match_pyrad_ann_df, 
                                              feat_names = feature_names, 
                                              epsilons = eps)
    
    annotation_pyrad_stats.to_csv(out_path / "recist_vs_diameter_stats.csv")


