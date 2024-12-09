# -*- coding: utf-8 -*-
"""RF_google_colab.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1rLA-3BDFB6XXJ28aV5XZmgDiHrQ0t-lW
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_validate
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
import joblib

# Load training data
df = pd.read_csv('QN_TCGA.csv', index_col=0)

df.columns

# Split the data into training and testing sets
X = df.iloc[:, :-1]
y = df['response']

X.shape

# Split the dataset into training (70%), test (20%), and validation (10%)
X_temp, X_test, y_temp, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.125, random_state=42, stratify=y_temp)  # 0.125 of 0.8 is 10%

X_train.shape

X_test.shape

X_val.shape

# 1. Standardization
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_val_scaled = scaler.transform(X_val)
X_test_scaled = scaler.transform(X_test)

smote = SMOTE(random_state=42)
rf_model = RandomForestClassifier(n_estimators=500,random_state=42)

pipeline = Pipeline(steps=[('smote', smote), ('rf_model', rf_model)])

skf =StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Perform cross-validation with CV = 10
cv_scores = cross_validate(pipeline, X_train_scaled, y_train, cv=skf, scoring ='accuracy', n_jobs=-1, return_estimator=True, return_train_score=True)

# Save the model to a file
joblib.dump(cv_scores, 'RF_New_QN_batch_TCGA_597.pkl')
from google.colab import files
files.download('RF_New_QN_batch_TCGA_597.pkl')

cv_scores

import matplotlib.pyplot as plt
import seaborn as sns
# Extract metrics
scores = cv_scores['test_score']

# Create a DataFrame for better handling
import pandas as pd
df = pd.DataFrame(scores, columns=['Accuracy'])

# Plot the scores
plt.figure(figsize=(10, 6))
sns.boxplot(data=df)
plt.title('Cross-Validation Scores')
plt.ylabel('Accuracy')
plt.show()

fit_times = cv_scores['fit_time']
score_times = cv_scores['score_time']

# Create DataFrames
df_fit = pd.DataFrame(fit_times, columns=['Fit Time'])
df_score = pd.DataFrame(score_times, columns=['Score Time'])

# Plot training time
plt.figure(figsize=(10, 6))
sns.boxplot(data=df_fit)
plt.title('Training Time Across Folds')
plt.ylabel('Time (seconds)')
plt.show()

# Plot scoring time
plt.figure(figsize=(10, 6))
sns.boxplot(data=df_score)
plt.title('Scoring Time Across Folds')
plt.ylabel('Time (seconds)')
plt.show()

# Print training and test scores
print("Training scores for each fold:", cv_scores['train_score'])
print("Test scores for each fold:", cv_scores['test_score'])

# Calculate and print the average cross-validation score
print("Average cross-validation score:", cv_scores['test_score'].mean())

cv_scores['estimator']

from sklearn.metrics import accuracy_score
for idx, estimator in enumerate(cv_scores['estimator']):
    print(f"\nModel {idx+1} trained on fold {idx+1}")
    # You can further inspect or use the estimator, e.g., predict on the test set
    y_pred_fold = estimator.predict(X_test_scaled)
    print(f"Accuracy on test set with model {idx+1}: {accuracy_score(y_test, y_pred_fold):.4f}")

pipeline.fit(X_train_scaled, y_train)

# Save the model to a file
joblib.dump(pipeline, 'RF_new_pipeline_QN_TCGA_batch_597.pkl')
#from google.colab import files
#files.download('rf_model_piepline_597.pkl')

y_pred = pipeline.predict(X_test_scaled)

y_pred_val = pipeline.predict(X_val_scaled)

# Assuming y_test and y_pred are already defined
results_df = pd.DataFrame({
    'y_test': y_test,
    'y_pred': y_pred
})

# Save to CSV
results_df.to_csv('predictions_RF_new_QN_TCGA_batch_597.csv', index=False)

# Assuming y_test and y_pred are already defined
results_df = pd.DataFrame({
    'y_val': y_val,
    'y_pred': y_pred_val
})

# Save to CSV
results_df.to_csv('predictions_RFQNTCGAbatch597_intVal.csv', index=False)

from sklearn.metrics import classification_report
report = classification_report(y_test, y_pred, zero_division=0, output_dict=True)
# Save classification report for boruta features to CSV
report_df = pd.DataFrame(report).transpose()
report_df.to_csv('report_RFnew_QN_TCGA_batch_597.csv', index=True)

report_val = classification_report(y_val, y_pred_val, zero_division=0, output_dict=True)
# Save classification report for boruta features to CSV
report_df_val = pd.DataFrame(report_val).transpose()
report_df_val.to_csv('report_RFQNTCGAbatch597_intVal.csv', index=True)

# prompt: give me confusion matrix plot and the lables take from y_test

import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

cm = confusion_matrix(y_test, y_pred)

plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=np.unique(y_test), yticklabels=np.unique(y_test))
plt.title('Confusion Matrix')
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.savefig('confusion_matrix_RF_QN_TCGA_batch_597.jpeg', dpi=300, format='jpeg')
plt.show()

cm_val = confusion_matrix(y_val, y_pred_val)

plt.figure(figsize=(8, 6))
sns.heatmap(cm_val, annot=True, fmt='d', cmap='Blues',
            xticklabels=np.unique(y_val), yticklabels=np.unique(y_val))
plt.title('Confusion Matrix')
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.savefig('cm_RFQNTCGAbatch597_intVal.jpeg', dpi=300, format='jpeg')
plt.show()

from sklearn.metrics import roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
# Confusion Matrix
conf_matrix = confusion_matrix(y_test, y_pred)
conf_matrix_df = pd.DataFrame(conf_matrix, index=np.unique(y_test), columns=np.unique(y_test))
conf_matrix_df.to_csv('cm_RF_TCGA_QN_batch_597.csv',float_format='%.0f')

# Confusion Matrix
conf_matrix_val = confusion_matrix(y_val, y_pred_val)
conf_matrix_df_val = pd.DataFrame(conf_matrix_val, index=np.unique(y_val), columns=np.unique(y_val))
conf_matrix_df_val.to_csv('cm_RFTCGAQNbatch597_intVal.csv',float_format='%.0f')

y_prob = pipeline.predict_proba(X_test_scaled)  # For ROC curve

y_prob_val = pipeline.predict_proba(X_val_scaled)

from sklearn.preprocessing import label_binarize
# Assuming you have a mapping from class indices to actual class names
class_names = np.unique(y_test)  # This gets the actual class names
# Binarize the test labels
y_test_binarized = label_binarize(y_test, classes=class_names)
joblib.dump(y_test_binarized, 'RF_TCGA_QN_batch_597_y_binary.pkl')

joblib.dump(y_prob, 'RF_TCGA_QN_batch_597_y_prob.pkl')

# Assuming you have a mapping from class indices to actual class names
class_names_val = np.unique(y_val)  # This gets the actual class names
# Binarize the test labels
y_test_bin_val = label_binarize(y_val, classes=class_names_val)
joblib.dump(y_test_bin_val, 'RFTCGAQNbatch597_intVal_ybinary.pkl')

joblib.dump(y_prob_val, 'RFTCGAQNbatch597_intVal_yprob.pkl')

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
plt.savefig('roc_RF_TCGA_QN_batch_597.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

# Plot legend in a separate figure
fig, ax = plt.subplots(figsize=(6, 6))

# Plot dummy lines to create a legend
for class_name in class_names:
    ax.plot([], [], lw=2, label=f'{class_name} (AUC = {roc_auc_dict[class_name]:0.2f})')

ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5), fontsize='large')
ax.axis('off')  # Hide axes

# Save the legend plot
plt.savefig('roc_legend_RF_TCGA_QN_batch_597.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

# Convert 180mm to inches (since Matplotlib uses inches)
#width_in_inches = 6  # 1 inch = 25.4 mm

#print("ROC and legend plots saved as 'roc_rf_rf.png' and 'roc_legend_rf_rf.png'.")

# Assuming you have a mapping from class indices to actual class names
class_names = np.unique(y_val)  # This gets the actual class names

# Binarize the test labels
y_test_bin_val = label_binarize(y_val, classes=class_names)


# Initialize lists to store ROC data
fpr_dict = {}
tpr_dict = {}
roc_auc_dict = {}

# Compute ROC curve and ROC area for each class
for i in range(y_test_bin_val.shape[1]):
    fpr, tpr, _ = roc_curve(y_test_bin_val[:, i], y_prob_val[:, i])
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
plt.savefig('rocRFTCGAQNbatch597_intVal.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

# Plot legend in a separate figure
fig, ax = plt.subplots(figsize=(6, 6))

# Plot dummy lines to create a legend
for class_name in class_names:
    ax.plot([], [], lw=2, label=f'{class_name} (AUC = {roc_auc_dict[class_name]:0.2f})')

ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5), fontsize='large')
ax.axis('off')  # Hide axes

# Save the legend plot
plt.savefig('roclegRFTCGAQNbatch597_intVal.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

# Convert 180mm to inches (since Matplotlib uses inches)
#width_in_inches = 6  # 1 inch = 25.4 mm

#print("ROC and legend plots saved as 'roc_rf_rf.png' and 'roc_legend_rf_rf.png'.")

# Assuming you have fitted your model and have predictions
train_accuracy = pipeline.score(X_train_scaled, y_train)
test_accuracy = pipeline.score(X_test_scaled, y_test)
val_accuracy = pipeline.score(X_val_scaled, y_val)

print("Training Accuracy:", train_accuracy)
print("Test Accuracy:", test_accuracy)
print("Validation Accuracy:", val_accuracy)

if train_accuracy - val_accuracy > 0.1:  # Example threshold
    print("Model may be overfitting.")
else:
    print("Model performance is consistent across train and validation sets.")

if train_accuracy - test_accuracy > 0.1:  # Example threshold
    print("Model may be overfitting.")
else:
    print("Model performance is consistent across train and test sets.")

import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import learning_curve

# Assuming you have your pipeline defined as:
# pipeline = Pipeline(steps=[('smote', smote), ('rf_model', rf_model)])

# Generate learning curves
train_sizes, train_scores, test_scores = learning_curve(
    pipeline,
    X_train_scaled,  # Use the original training data before SMOTE
    y_train,
    train_sizes=np.linspace(0.1, 1.0, 10),  # Use 10 points from 10% to 100% of the training data
    cv=skf,  # 5-fold cross-validation
    n_jobs=-1,  # Use all available cores
    scoring='accuracy'  # You can change this to another metric if needed
)

# Calculate mean and standard deviation for training scores
train_scores_mean = np.mean(train_scores, axis=1)
train_scores_std = np.std(train_scores, axis=1)

# Calculate mean and standard deviation for validation scores
test_scores_mean = np.mean(test_scores, axis=1)
test_scores_std = np.std(test_scores, axis=1)

# Plotting the learning curves
plt.figure(figsize=(10, 6))
plt.title('Learning Curves (Random Forest with SMOTE)')
plt.xlabel('Training Size')
plt.ylabel('Accuracy')
plt.ylim(0, 1)
plt.grid()

# Plot training scores
plt.plot(train_sizes, train_scores_mean, 'o-', color='r', label='Training Score')
plt.fill_between(train_sizes,
                 train_scores_mean - train_scores_std,
                 train_scores_mean + train_scores_std,
                 color='r', alpha=0.1)

# Plot validation scores
plt.plot(train_sizes, test_scores_mean, 'o-', color='g', label='Validation Score')
plt.fill_between(train_sizes,
                 test_scores_mean - test_scores_std,
                 test_scores_mean + test_scores_std,
                 color='g', alpha=0.1)

plt.legend(loc='best')
plt.show()



# Loop through each class to plot the ROC curve separately
for i, class_name in enumerate(class_names):
    # Create a new figure for each class
    plt.figure(figsize=(8, 6))

    # Plot ROC curve for the current class
    plt.plot(fpr_dict[class_name], tpr_dict[class_name], lw=2, color='blue', label=f'ROC curve (AUC = {roc_auc_dict[class_name]:0.2f})')

    # Plot the diagonal line for random guessing
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', lw=2)

    # Set the plot limits and labels
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve for {class_name}')

    # Add the legend
    plt.legend(loc='lower right')

    # Save the plot for the current class
    plt.savefig(f'roc_curve_{class_name}.png', bbox_inches='tight')
    plt.close()