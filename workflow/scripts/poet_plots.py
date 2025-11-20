from damply import dirs
import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np 
from pathlib import Path

np.random.seed(10)

# N = 100
# dice_scores = np.random.rand(N)
# hd95_scores = np.random.uniform(0, 250, N)

# recist_class_flip_bool = np.random.choice(a=[False, True], size=N, p=[0.3, 0.7])

dataset = "TCIA_CPTAC-PDA"

results_df = pd.read_csv(dirs.RESULTS / "POETPresentation" / f"indiv_{dataset}_results.csv", index_col=0)

out_dir = Path(dirs.RESULTS / "POETPresentation" / "plots" / dataset)
out_dir.mkdir(parents=True, exist_ok=True)

dice_scores = results_df['volume_dice']
hd95_scores = results_df['hausdorff']

recist_class_flip_bool = results_df['PredRECISTCat'] == results_df['RECISTCat']


# print(len(recist_class_flip_bool == True))
dice_same = dice_scores[recist_class_flip_bool]
hd95_same = hd95_scores[recist_class_flip_bool]
print(len(dice_same))

dice_diff = dice_scores[~recist_class_flip_bool]
hd95_diff = hd95_scores[~recist_class_flip_bool]
print(len(dice_diff))
fig, ax = plt.subplots()

scatter_same = ax.scatter(dice_same, hd95_same, c='blue', marker='o', label='Same RECIST class')
scatter_diff = ax.scatter(dice_diff, hd95_diff, c='red', marker='^', label='Changed RECIST class')
ax.set_xlabel('DICE score')
ax.set_ylabel('Hausdorff 95 distance')
fig.suptitle('MedSAM2-RECIST Prediction Evaluation')
ax.set_title(dataset)
ax.legend(loc='upper right')

plt.savefig(out_dir / "dice_v_hd95_recist_flip.png")


fig2, ax2 = plt.subplots()

vol_scatter = ax2.scatter(results_df['T0VoxVol'], results_df['predVoxVol'])
ax2.set_xlabel('Ground Truth Voxel Volume')
ax2.set_ylabel('Predicted Voxel Volume')
fig.suptitle('MedSAM2-RECIST Volume Comparison')
ax2.set_title(dataset)

plt.savefig(out_dir / "orig_vs_pred_voxvol.png")


fig3, ax3 = plt.subplots(1,8, layout='constrained')
fig3.set_size_inches(40,5)

metrics = ['volume_dice','jaccard','hausdorff','surface_dice','panoptic_quality','added_path_length','false_negative_volume','false_negative_path_length']
colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'orange', 'purple', 'black']

fig.suptitle(dataset)
for idx, metric in enumerate(metrics):
    ax3[idx].scatter(results_df['T0VoxVol'], results_df[metric], c=colors[idx])
    ax3[idx].set_title(metric)
    ax3[idx].set_ylabel(f"{metric} vs. Ground Truth Volume")
    ax3[idx].set_xlabel('Ground Truth Voxel Volume')

plt.savefig(out_dir / "metrics_v_volume.png")


    