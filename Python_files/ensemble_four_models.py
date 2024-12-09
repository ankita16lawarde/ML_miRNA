from sklearn.ensemble import VotingClassifier, RandomForestClassifier, AdaBoostClassifier
from xgboost import XGBClassifier
import lightgbm as lgb
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, accuracy_score,roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
from imblearn.pipeline import Pipeline  # Note the use of imbalanced-learn's pipeline
from imblearn.over_sampling import SMOTE
from sklearn.tree import DecisionTreeClassifier
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("/gpfs/helios/home/etais/hpc_ankita.ee/miRNA_network_ML/QN_TCGA.csv", index_col = 0)

X = df.iloc[:,:-1]
y = df['response']

#RFE_features = pd.read_csv('RFE_features_150new.csv')
#RF_features = pd.read_csv('rf_selected_features.csv')
#boruta = pd.read_csv('boruta_selected_features.csv')
#lda = pd.read_csv('lda_selected_features.csv')


#RFE_features = RFE_features['0']
#RF_features = RF_features['0']
#boruta_features = boruta['0']
#lda_features = lda['0']

#RFE_features = RFE_features.tolist()
#RF_features = RF_features.tolist()
#boruta_features = boruta_features.tolist()
#lda_features = lda_features.tolist()

#RFE_subset = X[X.columns.intersection(RFE_features)]
#rf_subset = X[X.columns.intersection(RF_features)]
#boruta_subset = X[X.columns.intersection(boruta_features)]
#lda_subset = X[X.columns.intersection(lda_features)]

# Split the dataset into training (70%), test (20%), and validation (10%)
X_temp, X_test, y_temp, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.125, random_state=42, stratify=y_temp)  # 0.125 of 0.8 is 10%

#X_temp, X_test, y_temp, y_test = train_test_split(lda_subset, y, test_size=0.2, random_state=42, stratify=y)
#X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.125, random_state=42, stratify=y_temp)  # 0.125 of 0.8 is 10%

# Assuming 'y' is a pandas Series or a numpy array and 'subset_LDA' is a pandas DataFrame or numpy array
class_to_remove = 'TCGA-ESCA-NT'

# Create a boolean mask for filtering
mask_train = y_train != class_to_remove
mask_val = y_val != class_to_remove
mask_test = y_test != class_to_remove

# Filter the datasets
X_train_filtered = X_train[mask_train]
y_train_filtered = y_train[mask_train]

X_val_filtered = X_val[mask_val]
y_val_filtered = y_val[mask_val]

X_test_filtered = X_test[mask_test]
y_test_filtered = y_test[mask_test]

# Now, X_train_filtered, y_train_filtered, X_val_filtered, y_val_filtered, X_test_filtered, y_test_filtered
# contain the data without the TCGA-ESCA-NT class

# Convert string labels to numerical labels
label_encoder = LabelEncoder()
y_train_filtered = label_encoder.fit_transform(y_train_filtered)
y_val_filtered = label_encoder.fit_transform(y_val_filtered)
y_test_filtered = label_encoder.fit_transform(y_test_filtered)

# Combine the filtered validation and test sets
X_combined = pd.concat([X_val_filtered, X_test_filtered], ignore_index=True)
y_combined = pd.concat([pd.Series(y_val_filtered), pd.Series(y_test_filtered)], ignore_index=True)

# Step 2: Initialize individual classifiers
rf = RandomForestClassifier(n_estimators=100, random_state=42)

xgboost = XGBClassifier(objective='multi:softmax', num_class=28, random_state=42)

ada = AdaBoostClassifier(base_estimator = DecisionTreeClassifier(max_depth=5),
                         n_estimators=300,
                         learning_rate=1.2,
                         algorithm="SAMME")

lgbm = lgb.LGBMClassifier(random_state=42)

# Step 3: Create a VotingClassifier with soft voting
ensemble = VotingClassifier(
    estimators=[('rf', rf), ('xgboost', xgboost), ('ada', ada), ('lgbm', lgbm)],
    voting='soft'
)

pipeline = Pipeline([
    ('scaler', StandardScaler()),  # Standardize features
    ('smote', SMOTE(random_state=42)),  # Handle class imbalance
    ('ensemble', ensemble)  # Xensemble of fours methods
])

# Initialize Stratified K-Fold
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Perform cross-validation
cv_results = cross_validate(pipeline, X_train_filtered, y_train_filtered, cv=skf, scoring='accuracy', return_train_score=True, n_jobs=-1, return_estimator=True)

#import joblib
# Save the model to a file
joblib.dump(cv_results, 'ensemble_model.pkl')
#from google.colab import files
#files.download('adaboost_597.pkl')

pipeline.fit(X_train_filtered, y_train_filtered)
# Save the model to a file
joblib.dump(pipeline, 'ensemble_pipeline.pkl')

y_pred = pipeline.predict(X_combined)
#y_pred_val = pipeline.predict(X_val)

# Assuming y_test and y_pred are already defined
results_df = pd.DataFrame({
    'y_test': y_combined,
    'y_pred': y_pred
})

# Save to CSV
results_df.to_csv('predictions_ensemble.csv', index=False)


report = classification_report(y_combined, y_pred, zero_division=0, output_dict=True)
# Save classification report for boruta features to CSV
report_df = pd.DataFrame(report).transpose()
report_df.to_csv('report_ensemble.csv', index=True)

#report_val = classification_report(y_val, y_pred_val, zero_division=0, output_dict=True)
# Save classification report for boruta features to CSV
#report_df_val = pd.DataFrame(report_val).transpose()
#report_df_val.to_csv('report_adaLDA_intVal.csv', index=True)


cm = confusion_matrix(y_combined, y_pred)

plt.figure(figsize=(10, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=np.unique(y_combined), yticklabels=np.unique(y_combined))
plt.title('Confusion Matrix')
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.tight_layout()  # Adjust layout
plt.savefig('cm_ensemble.jpeg', dpi=300, format='jpeg', bbox_inches='tight')
#plt.show()

#cm_val = confusion_matrix(y_val, y_pred_val)

#plt.figure(figsize=(10, 6))
#sns.heatmap(cm_val, annot=True, fmt='d', cmap='Blues',
  #          xticklabels=np.unique(y_val), yticklabels=np.unique(y_val))
#plt.title('Confusion Matrix')
#plt.xlabel('Predicted Label')
#plt.ylabel('True Label')
#plt.tight_layout()  # Adjust layout
#plt.savefig('cm_adaLDA_intVal.jpeg', dpi=300, format='jpeg')
#plt.show()


# Confusion Matrix
conf_matrix = confusion_matrix(y_combined, y_pred)
conf_matrix_df = pd.DataFrame(conf_matrix, index=np.unique(y_combined), columns=np.unique(y_combined))
conf_matrix_df.to_csv('cm_ensemble.csv',float_format='%.0f')

# Confusion Matrix
#conf_matrix_val = confusion_matrix(y_val, y_pred_val)
#conf_matrix_df_val = pd.DataFrame(conf_matrix_val, index=np.unique(y_val), columns=np.unique(y_val))
#conf_matrix_df_val.to_csv('cm_adaLDA_intVal.csv',float_format='%.0f')

y_prob = pipeline.predict_proba(X_combined)  # For ROC curve

#y_prob_val = pipeline.predict_proba(X_val)


# Assuming you have a mapping from class indices to actual class names
class_names = np.unique(y_combined)  # This gets the actual class names
# Binarize the test labels
y_test_binarized = label_binarize(y_combined, classes=class_names)
joblib.dump(y_test_binarized, 'ensemble_y_binary.pkl')

joblib.dump(y_prob, 'ensemble_y_prob.pkl')

# Assuming you have a mapping from class indices to actual class names
#class_names_val = np.unique(y_val)  # This gets the actual class names
# Binarize the test labels
#y_test_bin_val = label_binarize(y_val, classes=class_names_val)
#joblib.dump(y_test_bin_val, 'adaLDA_intVal_ybinary.pkl')

#joblib.dump(y_prob_val, 'adaLDA_intVal_yprob.pkl')


# Assuming you have a mapping from class indices to actual class names
class_names = np.unique(y_combined)  # This gets the actual class names

# Binarize the test labels
y_test_binarized = label_binarize(y_combined, classes=class_names)


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
plt.savefig('roc_ensemble.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()

# Plot legend in a separate figure
fig, ax = plt.subplots(figsize=(6, 6))

# Plot dummy lines to create a legend
for class_name in class_names:
    ax.plot([], [], lw=2, label=f'{class_name} (AUC = {roc_auc_dict[class_name]:0.2f})')

ax.legend(loc='center', bbox_to_anchor=(0.5, 0.5), fontsize='large')
ax.axis('off')  # Hide axes

# Save the legend plot
plt.savefig('roc_legend_ensemble.jpg', dpi=300, format='jpeg', bbox_inches='tight')
plt.close()