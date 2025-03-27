import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, accuracy_score
from imblearn.pipeline import Pipeline  # Note the use of imbalanced-learn's pipeline
from imblearn.over_sampling import SMOTE
import lightgbm as lgb

# Load the HDF5 file into a Pandas dataframe
#df = pd.read_csv("QN_TCGA.csv", index_col = 0)
df = pd.read_csv("protein_coding_matrix_RNA-Seq.csv", index_col= 0)

#df

# Split features and target
X = df.iloc[:, :-1]  # Features
y = df['response']  # Target variable


# Split the dataset into training (70%), test (20%), and validation (10%)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
#X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.125, random_state=42, stratify=y_temp)  # 0.125 of 0.8 is 10%

# Define the pipeline with StandardScaler, SMOTE, and LightGBM
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('smote', SMOTE(random_state=42)),
    ('lightgbm', lgb.LGBMClassifier(random_state=42))
])

# Initialize Stratified K-Fold
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Perform cross-validation
cv_results = cross_validate(pipeline, X_train, y_train, cv=skf, scoring='accuracy', return_train_score=True, n_jobs=-1, return_estimator=True)

import joblib
joblib.dump(cv_results, 'lightGBM_model.pkl')

# Print the cross-validation results
print("Training Accuracy (mean):", np.mean(cv_results['train_score']))
print("Validation Accuracy (mean):", np.mean(cv_results['test_score']))


# Fit the model on the full training data
pipeline.fit(X_train, y_train)

joblib.dump(pipeline, 'lightGBM_pipeline.pkl')


y_pred = pipeline.predict(X_test)

# Assuming y_test and y_pred are already defined
results_df = pd.DataFrame({
    'y_test': y_test,
    'y_pred': y_pred
})

# Save to CSV
results_df.to_csv('predictions_lightGBM.csv', index=False)

from sklearn.metrics import classification_report
report = classification_report(y_test, y_pred, output_dict=True)

# Save classification report for boruta features to CSV
report_df = pd.DataFrame(report).transpose()
report_df.to_csv('classification_report_lightGBM.csv', index=True)


# prompt: give me confusion matrix plot and the lables take from y_test
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

cm = confusion_matrix(y_test, y_pred)

plt.figure(figsize=(10, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=np.unique(y_test), yticklabels=np.unique(y_test))
plt.title('Confusion Matrix')
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.savefig('cm_lightGBM.jpeg', dpi=300, format='jpeg', bbox_inches='tight')
plt.show()

from sklearn.metrics import roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
# Confusion Matrix
conf_matrix = confusion_matrix(y_test, y_pred)
conf_matrix_df = pd.DataFrame(conf_matrix, index=np.unique(y_test), columns=np.unique(y_test))
conf_matrix_df.to_csv('cm_lightGBM.csv',float_format='%.0f')

y_prob = pipeline.predict_proba(X_test)  # For ROC curve

#y_prob_val = pipeline.predict_proba(X_val)  # For ROC curve

from sklearn.preprocessing import label_binarize
# Assuming you have a mapping from class indices to actual class names
class_names = np.unique(y_test)  # This gets the actual class names
# Binarize the test labels
y_test_binarized = label_binarize(y_test, classes=class_names)
joblib.dump(y_test_binarized, 'lightGBM_y_binary.pkl')

joblib.dump(y_prob, 'lightGBM_y_prob.pkl')

import matplotlib.pyplot as plt

# Assuming you have a mapping from class indices to actual class names
class_names = np.unique(y_test)  # This gets the actual class names

# Binarize the test labels
y_test_binarized = label_binarize(y_test, classes=class_names)


# Initialize lists to store ROC data
fpr_dict = {}
tpr_dict = {}
roc_auc_dict = {}

# Compute ROC curve and ROC area for each class
for i in range(y_test_binarized.shape[1]):
    fpr, tpr, _ = roc_curve(y_test_binarized[:, i], y_prob[:, i])
    roc_auc = auc(fpr, tpr)
    fpr_dict[class_names[i]] = fpr
    tpr_dict[class_names[i]] = tpr
    roc_auc_dict[class_names[i]] = roc_auc

# Plot ROC curves
plt.figure(figsize=(10, 10))
for class_name in class_names:
    plt.plot(fpr_dict[class_name], tpr_dict[class_name], lw=2, label=f'{class_name} (AUC = {roc_auc_dict[class_name]:0.2f})')

plt.plot([0, 1], [0, 1], color='gray', linestyle='--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve for Multiclass Classification')

# Save the ROC plot without the legend
plt.savefig('roc_lightGBM.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

# Plot legend in a separate figure
fig, ax = plt.subplots(figsize=(6, 6))

# Plot dummy lines to create a legend
for class_name in class_names:
    ax.plot([], [], lw=2, label=f'{class_name} (AUC = {roc_auc_dict[class_name]:0.2f})')

ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5), fontsize='large')
ax.axis('off')  # Hide axes

# Save the legend plot
plt.savefig('roc_legend_lightGBM.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

