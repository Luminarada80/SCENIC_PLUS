import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc
import matplotlib.pyplot as plt

# Step 1: Load the TSV files
inferred_df = pd.read_csv('/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/updated_inferred_network.tsv', sep='\t', low_memory=False)
ground_truth_df = pd.read_csv('/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/data/Ground_truth_mESC_RN111.tsv', sep='\t', quoting=3)
true_positives_df = pd.read_csv('/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/merged_data.csv')

# Step 2: Prepare true labels and predicted scores
# Create a set of true edges from the ground truth
true_edges = set(zip(ground_truth_df['Source'], ground_truth_df['Target']))

# Convert importance_score to numeric, forcing errors to NaN and then dropping NaNs
inferred_df['importance'] = pd.to_numeric(inferred_df['importance'], errors='coerce')
inferred_df.dropna(subset=['importance'], inplace=True)

# Create a DataFrame for the inferred network
inferred_edges = set(zip(inferred_df['TF'], inferred_df['target']))

# Create true labels and scores
y_true = []
y_scores = []

for index, row in inferred_df.iterrows():
    edge = (row['TF'], row['target'])  # Adjust these based on your actual column names
    if edge in true_edges:
        y_true.append(1)  # True positive
        y_scores.append(row['importance'])  # Score from inferred network
    else:
        y_true.append(0)  # False positive
        y_scores.append(row['importance'])  # Score from inferred network

# Ensure arrays are numpy arrays
y_true = np.array(y_true)
y_scores = np.array(y_scores)

# Step 3: Calculate AUROC
auroc = roc_auc_score(y_true, y_scores)
print(f"AUROC: {auroc:.2f}")

# Step 4: Plot the ROC curve
fpr, tpr, thresholds = roc_curve(y_true, y_scores)

plt.figure()
plt.plot(fpr, tpr, color='blue', label='ROC curve (area = %0.2f)' % auroc)
plt.plot([0, 1], [0, 1], color='red', linestyle='--')  # Diagonal line
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc='lower right')
plt.show()

# Step 5: Calculate AUPRC
precision, recall, _ = precision_recall_curve(y_true, y_scores) 
auprc = auc(recall, precision)
print(f"AUPRC: {auprc:.2f}")

# Step 6: Plot the Precision-Recall curve
plt.subplot(1, 2, 2) # Precision-Recall Curve
plt.plot(recall, precision, color='green', label='Precision-Recall curve (area = %0.2f)' % auprc)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc='lower left')
plt.tight_layout()
plt.show()
