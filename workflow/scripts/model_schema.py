### Prototype hybrid class for AAuRA segmentation models ###
import numpy as np
import pandas as pd

from pathlib import Path
from utils.utils import calc_metrics

class SegModel(): 
    '''  
    Hybrid class for the initialization, inference, and analysis/post-processing of outputs for various segmentation
    models. 
    '''
    def __init__(self, 
                 model_type: str):
        self.model_type = model_type
        
        ## Necessary properties for methods below ##
        #To be definitively defined through the created subclasses as the 
        #pipeline is ran
        # Paths
        self.img_path: Path = None # Imaging path 
        self.gts_path: Path = None # Ground truth segmentation path
        self.checkpoint_path: Path = None #Where the model checkpoint is stored (to be used for init_model)
        self.prompt_skel_path: Path = None # For BiomedParse text prompt skeletons

        # Imaging and Segmentation Arrays
        self.img_arr: np.ndarray = None # Imaging array
        self.gts_arr: np.ndarray = None # Ground truth segmentation array (binary)
        self.pred_arr: np.ndarray = None # Predicted segmentation array 

        # Image Properties 
        self.spacing: np.ndarray = None 
        self.orientation: np.ndarray = None
        self.direction: np.ndarray = None

        # Preprocessing Properties 
        self.preproc_info: dict = None # What preprocessing was used (to be used for visualizations)

    ### Abstract methods ###
    def init_model(self, 
                   checkpoint_path: Path):
        ''' 
        Used to pass in model checkpoint. Must be implemented by 
        subclass.

        Parameters
        ----------
        checkpoint_path: Path
            Location of the model checkpoint 
        '''
        raise NotImplementedError()
        
    def inference(self): 
        '''  
        Perform auto-segmentation on image based on prompt. Must be implemented
        by subclass. 
        To be run after img_arr, gts_arr, and necessary prompt information 
        has been added to the object instance.
        '''
        raise NotImplementedError() 
    
    ### Concrete methods ###
    def get_metrics(self) -> pd.DataFrame:
        '''
        Calculate performance metrics based on the predicted and ground truth masks and save into a dataframe.
        '''
        metric_df = calc_metrics(pred_mask = self.pred_arr, 
                                 gt_mask = self.gts_arr, 
                                 spacing = self.spacing, 
                                 filename = self.gts_path,
                                 text_prompt = self.text_prompt)
        return metric_df
    
    def visualizations(self): 
        '''
        
        '''