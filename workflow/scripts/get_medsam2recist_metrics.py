import numpy as np 
import pandas as pd

from pathlib import Path
from evaluate import Evaluator 
from damply import dirs

def compile_metrics_eff(results_path: Path, 
                    save_path: Path, 
                    ): 
    '''
    Calculate performance metrics of all generated segmentations from MedSAM2-RECIST models. Results must
    be in .npz format with the following keys: segs, gts, boxes, spacing

    Parameters
    ----------
    results_path: Path 
        Path to the .npz files that were created by MedSAM2-RECIST
    save_path: Path
        Where to save the csv with all metric results 
    '''
    # Loop over all of the data within specified folder 
    for file in results_path.iterdir(): 
        # Ignore files that are not npz 
        if not str(file).endswith('.npz'): 
            continue 
        
        # Load in npz data 
        seg_model_result = np.load(file)

        # Get data required for metrics calculation 
        gen_seg = seg_model_result['segs']
        grd_truth_seg = seg_model_result['gts']
        seg_spacing = seg_model_result['spacing']

        # Create instance of metric evaluator and get all metrics in the form of a dictionary 
        metric_eval = Evaluator()
        curr_metric_dict = metric_eval(preds = gen_seg, 
                                       targets = grd_truth_seg, 
                                       spacing = seg_spacing)
        
        # Convert dictionary to dataframe, add filename info, and add performance to overall dataframe 
        curr_metric_df = pd.DataFrame(curr_metric_dict, index = [0]) 
        
        # Save relative path starting from dataset folder 
        rel_file = '/'.join(str(file).split('/')[-4:])
        curr_metric_df['filename'] = rel_file

        if 'all_metrics_df' not in locals(): # Check if dataframe variable doesn't exist yet
            all_metrics_df = curr_metric_df
        else: 
            all_metrics_df = pd.concat([all_metrics_df, curr_metric_df], ignore_index = True).reset_index(drop = True) 
        
    # Check if save path exists
    if not save_path.exists(): 
        Path(save_path.mkdir(parents = True, exist_ok = True))

    save_name = save_path / "medsam2recisteff_seg_eval.csv"

    all_metrics_df.to_csv(save_name)

if __name__ == '__main__': 
    seg_path = Path('/home/bhkuser/bhklab/radiomics/Projects/MedSAM2-RECIST/data/results/RADCURE_OCSCC/predictions/eff-tiny')
    save_path = dirs.RESULTS / Path('RADCURE_OCSCC/metric_eval')
    compile_metrics_eff(results_path = seg_path, 
                    save_path = save_path)
        