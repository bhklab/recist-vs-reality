import numpy as np 
import pandas as pd

from pathlib import Path
from evaluate import Evaluator 
from damply import dirs

def compile_metrics(dataset: str, 
                    disease_site: str, 
                    save_path: Path): 
    '''
    Calculate performance metrics of all generated segmentations from MedSAM2-RECIST models. Results must
    be in .npz format with the following keys: segs, gts, boxes, spacing. If those are not in found in the results,
    it will search the for the matching input path and get the necessary missing data from there. 

    Parameters
    ----------
    dataset: str 
        The full name of the dataset that was processed and predicted on (e.g. TCIA_CPTAC-CCRCC)
    disease_site: str
        Where the disease is located anatomically (e.g. Abdomen).  
    save_path: Path
        Where to save the final dataframe containing all prediction evaluations for a given dataset
    '''
    results_path = dirs.RESULTS / Path(dataset) / Path('predictions')

    # Loop over all of the data within specified folder 
    for folder in results_path.iterdir(): 
        curr_model_path = results_path / folder
        for file in curr_model_path.iterdir(): 
            # Ignore files that are not npz 
            if not str(file).endswith('.npz'): 
                continue 
            
            # Load in npz data 
            seg_model_result = np.load(file)

            # Get data required for metrics calculation 
            gen_seg = seg_model_result['segs']
            try: 
                grd_truth_seg = seg_model_result['gts'] #Some models (like the non-efficient MedSAM2-RECIST) may not have this so need a check
            except KeyError: 
                curr_filename = str(file).split('/')[-1]
                input_file = dirs.PROCDATA / Path(disease_site) / Path(dataset) / Path('images/npz_' + dataset) / Path(curr_filename) #Get data from input if not found
                input_data = np.load(input_file)
                grd_truth_seg = input_data['gts']
                del input_data #Free space
            seg_spacing = seg_model_result['spacing']

            # Free space
            del seg_model_result

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

            # Save the model type 
            model_type = str(folder).split('/')[-1]
            curr_metric_df['model_name'] = model_type

            if 'all_metrics_df' not in locals(): # Check if dataframe variable doesn't exist yet
                all_metrics_df = curr_metric_df
            else: 
                all_metrics_df = pd.concat([all_metrics_df, curr_metric_df], ignore_index = True).reset_index(drop = True) 
            
            # Check if save path exists
            if not save_path.exists(): 
                save_path.mkdir(parents = True, exist_ok = True)

            save_name = save_path / "medsam2recist_seg_eval.csv"

            all_metrics_df.to_csv(save_name)

if __name__ == '__main__': 
    dataset = 'TCIA_CPTAC-CCRCC'
    disease_site = 'Abdomen'
    save_path = dirs.RESULTS / Path('TCIA_CPTAC-CCRCC/metric_eval')
    compile_metrics(dataset = dataset, 
                    disease_site = disease_site, 
                    save_path = save_path)
    
        