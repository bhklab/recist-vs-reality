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
def get_filename_hn1(filename): 
    new = filename.split("/")[-1]
    return new

#### RADCURE RESULTS #### 
out_dir = Path(dirs.RESULTS / "POETPresentation" / "plots")
out_dir.mkdir(parents=True, exist_ok=True)

# radcure_df = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/medsam2recist_seg_eval_radcure.csv')

# fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
# fig.suptitle('RADCURE RERECIST Performance (n = 2996)')
# ax1.violinplot(dataset = radcure_df['volume_dice'])
# ax1.set_ylabel('Volumetric DICE') 
# ax1.tick_params(labelbottom = False)
# ax1.set_xticks([])
# ax2.violinplot(dataset = radcure_df['hausdorff'])
# ax2.set_ylabel('Hausdorff Distance')
# ax2.tick_params(labelbottom = False)    
# ax2.set_xticks([])
# ax3.violinplot(dataset = radcure_df['added_path_length'])
# ax3.set_ylabel('Added Path Length')
# ax3.tick_params(labelbottom = False)
# ax3.set_xticks([])

# plt.tight_layout()
# plt.savefig(out_dir / "test_binned_RECIST_prop.png")

# ## RADCURE PLOT DICE VS VOLUME WITH DISEASE SITE 
# clinical_radcure = pd.read_excel('/home/bhkuser/bhklab/radiomics/PublicDatasets/srcdata/HeadNeck/TCIA_RADCURE/clinical/RADCURE_Clinical_v04_20241219.xlsx')
# volume = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/gts_vols.csv')
# radcure_metrics = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/medsam2recist_seg_eval_radcure.csv')

# radcure_metrics = radcure_metrics.rename(columns = {'filename': 'filenameLONG'})
# radcure_metrics['filename'] = radcure_metrics['filenameLONG'].apply(get_filename_hn1)



# ##### RECIST PLOTS ##### 
# ccrcc_metrics = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/medsam2recist_seg_eval_ccrcc_recist.csv')
# hn1_metrics = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/medsam2recist_seg_eval_hn1.csv')
# nsclc_metrics = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/medsam2recist_seg_eval_nsclcradgen_recist.csv')
# edge_colors = ['red', 'blue', 'purple']
# colors = ['lightcoral', 'lightblue', 'lightpink']
# vp_recist = plt.violinplot(dataset = [ccrcc_metrics['volume_dice'],
#                          hn1_metrics['volume_dice'], 
#                          nsclc_metrics['volume_dice']], 
#                          showmedians = True)
# for i, body in enumerate(vp_recist['bodies']): 
#     body.set_facecolor(colors[i])
#     body.set_edgecolor(edge_colors[i])

# plt.ylabel('Volumetric Dice')
# plt.xticks([1, 2, 3], ['CCRCC\n(n = 32)', 'CPTAC\n(n = 22)', 'NSCLC\n(n = 73)'])
# plt.tight_layout()
# plt.savefig(out_dir / 'recist_medsam2_dice_violin.png')

# ##### HN1 RERECIST PLOTTING ###
hn1_vols = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/gts_vols_hn1.csv')
hn1_evals = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/medsam2recist_seg_eval_hn1.csv')

hn1_evals = hn1_evals.rename(columns = {'filename': 'filenameLONG'})
hn1_evals['filename'] = hn1_evals['filenameLONG'].apply(get_filename_hn1)

hn1_total = pd.merge(hn1_vols, hn1_evals, on = 'filename', how = 'left')

fig, ax = plt.subplots()
fig.suptitle('Volumetric DICE Relationship with Ground Truth Volume\n HN1 (n = 137)')
ax.scatter(hn1_total['T0VoxVol'], hn1_total['volume_dice'])
ax.set_xlabel('Volume')
ax.set_ylabel('Volumetric DICE')
plt.tight_layout()

plt.savefig(out_dir / 'hn1_volume_vs_dice.png')

# fig, ax = plt.subplots()
# ax.violinplot(dataset = hn1_total['volume_dice'])
# ax.set_ylabel('Volumetric DICE') 
# ax.tick_params(labelbottom = False)
# ax.set_xticks([])
# fig.set_figwidth(3)
# plt.tight_layout()
# plt.savefig(out_dir / 'hn1_dice_violin.png')

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.suptitle('HN1 Performance (n = 137)')
ax1.violinplot(dataset = hn1_total['volume_dice'])
ax1.set_ylabel('Volumetric DICE') 
ax1.tick_params(labelbottom = False)
ax1.set_xticks([])
ax2.violinplot(dataset = hn1_total['hausdorff'])
ax2.set_ylabel('Hausdorff Distance')
ax2.tick_params(labelbottom = False)    
ax2.set_xticks([])
ax3.violinplot(dataset = hn1_total['added_path_length'])
ax3.set_ylabel('Added Path Length')
ax3.tick_params(labelbottom = False)
ax3.set_xticks([])

plt.tight_layout()
plt.savefig(out_dir / "test_hn1_violin_all.png")

##### SIM DATA #####

# dataset = "TCIA_HEAD-NECK-RADIOMICS-HN1"

# results_df = pd.read_csv(dirs.RESULTS / Path("POETPresentation") / Path(f"indiv_{dataset}_results.csv"), index_col=0)

# out_dir = Path(dirs.RESULTS / "POETPresentation" / "plots" / dataset)
# out_dir.mkdir(parents=True, exist_ok=True)

# dice_scores = results_df['volume_dice']
# hd95_scores = results_df['hausdorff']

# recist_class_flip_bool = results_df['PredRECISTCat'] == results_df['RECISTCat']

# # print(len(recist_class_flip_bool == True))
# dice_same = dice_scores[recist_class_flip_bool]
# hd95_same = hd95_scores[recist_class_flip_bool]
# print(len(dice_same))

# dice_diff = dice_scores[~recist_class_flip_bool]
# hd95_diff = hd95_scores[~recist_class_flip_bool]
# print(len(dice_diff))
# fig, ax = plt.subplots()
# scatter_same = ax.scatter(dice_same, hd95_same, c='blue', marker='o', label='Same RECIST class')
# scatter_diff = ax.scatter(dice_diff, hd95_diff, c='red', marker='^', label='Changed RECIST class')
# ax.set_xlabel('DICE score')
# ax.set_ylabel('Hausdorff 95 distance')
# fig.suptitle('MedSAM2-RECIST Prediction Evaluation')
# ax.set_title(dataset)
# ax.legend(loc='upper right')

# plt.savefig(out_dir / "dice_v_hd95_recist_flip.png")

##### GET TOTAL SIM RESULTS #####
out_dir = Path(dirs.RESULTS / "POETPresentation" / "plots" / "AllData")
out_dir.mkdir(parents=True, exist_ok=True)

hn1 = pd.read_csv('recist-vs-reality/data/results/POETPresentation/indiv_TCIA_HEAD-NECK-RADIOMICS-HN1_results.csv')
ccrcc = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/data/results/POETPresentation/indiv_TCIA_CPTAC-CCRCC_results.csv')
pda = pd.read_csv('/home/bhkuser/bhklab/kaitlyn/recist-vs-reality/data/results/POETPresentation/indiv_TCIA_CPTAC-PDA_results.csv')
radgen = pd.read_csv('recist-vs-reality/data/results/POETPresentation/indiv_TCIA_NSCLC-Radiogenomics_results.csv')
rad = pd.read_csv('recist-vs-reality/data/results/POETPresentation/indiv_TCIA_NSCLC-Radiomics_results.csv')

dataset_list = [hn1, ccrcc, pda, radgen, rad]
all_dataset = pd.concat(dataset_list)

dice_scores = all_dataset['volume_dice']
hd95_scores = all_dataset['hausdorff']

recist_class_flip_bool = all_dataset['PredRECISTCat'] == all_dataset['RECISTCat']

# print(len(recist_class_flip_bool == True))
dice_same = dice_scores[recist_class_flip_bool]
hd95_same = hd95_scores[recist_class_flip_bool]
print(len(dice_same))

dice_diff = dice_scores[~recist_class_flip_bool]
hd95_diff = hd95_scores[~recist_class_flip_bool]
print(len(dice_diff))
fig, ax = plt.subplots()
scatter_same = ax.scatter(dice_same, hd95_same, c='blue', marker='o', label='Same RECIST class', alpha = 0.7)
scatter_diff = ax.scatter(dice_diff, hd95_diff, c='red', marker='^', label='Changed RECIST class', alpha = 0.5)
ax.set_xlabel('DICE Score')
ax.set_ylabel('Hausdorff 95 Distance')
fig.suptitle('MedSAM2-RECIST Simulated Prediction Evaluation')
ax.set_title('(n = 824)')
ax.legend(loc='upper right')

plt.savefig(out_dir / "dice_v_hd95_recist_flip_all_data.png")

### Stacked bar chart ### 

# Making misclassification proportions vs DICE plot 
all_dataset['same_RECIST'] = all_dataset['PredRECISTCat'] == all_dataset['RECISTCat']
bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
labels = ['(0.4, 0.5]', '(0.5, 0.6]', '(0.9, 1.0]', '(0.6, 0.7]', '(0.7, 0.8]', '(0.8, 0.9]', '(0.3, 0.4]', '(0.2, 0.3]', '(0.1, 0.2]']
all_dataset['binned_RECIST'] = pd.cut(all_dataset['volume_dice'], bins = bins)
all_dataset['binned_RECIST'] = all_dataset['binned_RECIST'].astype(str)
all_dataset['binned_RECIST'] = all_dataset['binned_RECIST'].replace('nan', '(0.0, 0.1]')
test = all_dataset.sort_values(by = 'binned_RECIST') 
misclass_num = []
correct_num = []
for bin_cat in test['binned_RECIST'].astype(str).unique(): 
    #Subset dataframe 
    curr_bin = test[test['binned_RECIST'] == bin_cat]
    total_count = curr_bin.shape[0] 
    misclassified = curr_bin[curr_bin['same_RECIST'] == False]
    misclass_count = misclassified.shape[0] 
    correct = curr_bin[curr_bin['same_RECIST'] == True].shape[0]
    misclass_num.append(misclass_count) 
    correct_num.append(correct)

# plt.hist([all_dataset['binned_RECIST'], all_dataset['same_RECIST']], bins = 10, stacked = True, color = ['r', 'b'])

# plt.savefit('hist_all_data_misclassify.png')
fig, ax = plt.subplots()
fig.suptitle('Proportions of Misclassified RECIST Categories\n(n = 824)')
x = ['(0.0, 0.1]', '(0.1, 0.2]', '(0.2, 0.3]', '(0.3, 0.4]', '(0.4, 0.5]', '(0.5, 0.6]', '(0.6, 0.7]', '(0.7, 0.8]', '(0.8, 0.9]', '(0.9, 1.0]']
ax.bar(x, misclass_num, color='r')
ax.bar(x, correct_num, bottom = misclass_num, color = 'b')
ax.set_xlabel('Binned DICE Score') 
ax.set_ylabel('Number of Tumours Classified')
ax.tick_params(axis='x', labelrotation=45)
ax.legend(["False", "True"], title = 'Correctly Classified')
plt.tight_layout()
plt.savefig(out_dir / "test_binned_RECIST_prop_all.png")

#Stacked bar chart 

# scatter_same = ax.scatter(dice_same, hd95_same, c='blue', marker='o', label='Same RECIST class')
# scatter_diff = ax.scatter(dice_diff, hd95_diff, c='red', marker='^', label='Changed RECIST class')
# ax.set_xlabel('DICE score')
# ax.set_ylabel('Hausdorff 95 distance')

# fig.suptitle('MedSAM2-RECIST Prediction Evaluation')
# ax.set_title(dataset)
# ax.legend(loc='upper right')

# plt.savefig(out_dir / "dice_v_hd95_recist_flip.png")


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
    ax3[idx].scatter(results_df[metric], results_df['T0VoxVol'], c=colors[idx] )
    ax3[idx].set_title(f"{metric} vs. Ground Truth Volume")
    ax3[idx].set_xlabel(metric)
    ax3[idx].set_ylabel('Ground Truth Voxel Volume')

plt.savefig(out_dir / "metrics_v_volume.png")

# Making misclassification proportions vs DICE plot 
results_df['same_RECIST'] = results_df['PredRECISTCat'] == results_df['RECISTCat']
fig.suptitle(dataset) 
bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
results_df['binned_RECIST'] = pd.cut(results_df['volume_dice'], bins = bins)
misclass_prop_list = []
for bin_cat in results_df['binned_RECIST'].unique(): 
    #Subset dataframe 
    curr_bin = results_df[results_df['binned_RECIST'] == bin_cat]
    total_count = curr_bin.shape[0] 
    misclassified = curr_bin[curr_bin['same_RECIST'] == False]
    misclass_count = misclassified.shape[0] 
    misclass_prop = misclass_count / total_count
    misclass_prop_list.append(misclass_prop)


fig.suptitle(dataset)
plt.bar(results_df['binned_RECIST'].unique().to_list(), misclass_prop_list)
