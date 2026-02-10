import numpy as np
import pandas as pd

from damply import dirs
from evaluate import Evaluator
from pathlib import Path
from tqdm import tqdm


def run_evaluator(gt: np.ndarray, 
                  pred: np.ndarray,
                  spacing: tuple[float, float, float]
                  ) -> dict:
    """Run evaluation metrics on the predicted segementation mask against the ground truth mask.
    Args:
        gt_mask (np.ndarray): Ground truth segmentation mask.
        pred_mask (np.ndarray): Predicted segmentation mask.
    Returns:
        dict: Dictionary containing evaluation metric results."""
    evaluator = Evaluator(metrics = [
                                     "volume_dice", 
                                     "jaccard", 
                                     "hausdorff", 
                                     "surface_dice", 
                                     "panoptic_quality",
                                     "added_path_length", 
                                     "false_negative_volume", 
                                     "false_negative_path_length"
                                    ]
                         )
    
    evaluation_results = evaluator(preds = (pred == 1).astype(np.uint8),
                                   targets = (gt == 1).astype(np.uint8),
                                   spacing = spacing
                                   )
    return evaluation_results



def load_gt_npz(file_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load ground truth .npz file from MedSAM2-RECIST inference inputs.
    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Arrays for image, ground truth mask, recist annotation, spacing, origin, direction.
    """
    gt_data = np.load(file_path)
    image = gt_data["imgs"]
    gt = gt_data["gts"]
    recist = gt_data["recist"]
    spacing = gt_data["spacing"]
    origin = gt_data["origin"]
    direction = gt_data["direction"]
    return image, gt, recist, spacing, origin, direction



def load_pred_npz(file_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load prediction .npz file from MedSAM2-RECIST inference.
    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Arrays for predicted mask, ground truth mask, bounding boxes, and spacing.
    """
    pred_data = np.load(file_path)
    pred = pred_data["segs"]
    gt = pred_data["gts"]
    boxes = pred_data["boxes"]
    spacing = pred_data["spacing"]
    return pred, gt, boxes, spacing


if __name__ == "__main__":
    dataset = "RADCURE_OCSCC"

    sample_id_list = pd.read_csv(dirs.PROCDATA / dataset / "metadata" / "sample_ids.csv")["SampleID"].tolist()

    predictions_path = dirs.PROCDATA / dataset / "predictions" / "MedSAM2-RECIST" / "eff-tiny"

    results_df = []
    for sample_id in tqdm(sample_id_list, total=len(sample_id_list), desc="Evaluating samples"):
        print(f"Evaluating sample: {sample_id}")
        image, gt, recist, spacing, origin, direction = load_gt_npz(dirs.PROCDATA / dataset / "images" / f"npz_{dataset}" / f"{sample_id}.npz")
        pred, gt, boxes, spacing = load_pred_npz(predictions_path / f"{sample_id}.npz")

        sample_results = run_evaluator(gt=gt, pred=pred, spacing=tuple(spacing))

        results_df.append({"SampleID": sample_id,
                           **sample_results
                           })
        
        del image, gt, recist, spacing, origin, direction, pred, boxes, sample_results

    results_df = pd.DataFrame(results_df)
    evaluation_results_path = predictions_path / "sample_eval_metrics.csv"

    # average_results = 

    
