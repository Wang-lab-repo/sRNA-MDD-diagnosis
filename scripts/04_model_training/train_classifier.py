"""
train_classifier.py
===================
Random Forest classifier training with feature selection.

Pipeline:
  1. Load ZMD_GZ.csv from the data integration step
  2. Train full Random Forest for feature importance (gain + SHAP)
  3. Union top-15 from gain, SHAP (MDD vs others), and SHAP (other psychiatric diseases)
  4. Hyperparameter tuning via RandomizedSearchCV (5-fold stratified CV, AUC scoring)
  5. Train final model on all ZMD+GZ data with best parameters
  6. Save final_model.pkl and model_features.csv

Input:
  ./data/ZMD_GZ.csv

Output:
  ./result/train_model/final_model.pkl
  ./result/train_model/model_features.csv
  ./result/train_model/gain_importance.csv
  ./result/train_model/shap_importance.csv
  ./result/train_model/selected_features.csv
"""
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
import shap
import os
import re
import warnings
import joblib

warnings.filterwarnings("ignore")

# ================== Column name cleaner ==================
def clean_column_name(name):
    return re.sub(r'[^\w]', '_', str(name))

# ================== Paths ==================
data_path = "./data/ZMD_GZ.csv"          # from the preprocessing step
output_dir = "./result/train_model"
os.makedirs(output_dir, exist_ok=True)

# ================== 1. Load data ==================
df = pd.read_csv(data_path, index_col=0)
df.rename(columns=clean_column_name, inplace=True)

if "Group" not in df.columns:
    raise KeyError("Group column missing")
df["Group"] = pd.to_numeric(df["Group"], errors="coerce")
df = df.dropna(subset=["Group"])
df["Group"] = df["Group"].astype(int)

# Features start with miRNA_, rsRNA_, tsRNA_
feature_cols = [col for col in df.columns if col.startswith(('miRNA_', 'rsRNA_', 'tsRNA_'))]
print(f"Total features: {len(feature_cols)}")

# Label: MDD (Group==1) vs others
y = (df["Group"] == 1).astype(int)
X = df[feature_cols]

# ================== 2. Feature importance (full RF) ==================
print("Training full RF for feature importance...")
rf_full = RandomForestClassifier(
    n_estimators=200, random_state=42, class_weight="balanced", n_jobs=10
)
rf_full.fit(X, y)

# Gain importance
gain_imp = pd.DataFrame({
    "feature": feature_cols,
    "gain": rf_full.feature_importances_
}).sort_values("gain", ascending=False)

# SHAP importance
explainer = shap.TreeExplainer(rf_full)
shap_vals = explainer.shap_values(X)
if isinstance(shap_vals, list):
    shap_vals = shap_vals[1]
mean_abs_shap = np.abs(shap_vals).mean(axis=0)
shap_imp = pd.DataFrame({
    "feature": feature_cols,
    "shap": mean_abs_shap
}).sort_values("shap", ascending=False)

top15_gain = gain_imp["feature"].head(15).tolist()
top15_shap = shap_imp["feature"].head(15).tolist()
selected = sorted(set(top15_gain + top15_shap))
print(f"Initial features (union top15 gain/shap): {len(selected)}")

# ================== 3. Supplemental SHAP from other psychiatric diseases ==================
print("Computing SHAP for Group 2/3...")
mask_other = df["Group"].isin([2, 3])
if mask_other.sum() > 0:
    X_other = X[mask_other]
    shap_other = explainer.shap_values(X_other)
    if isinstance(shap_other, list):
        shap_other = shap_other[1]
    mean_abs_other = np.abs(shap_other).mean(axis=0)
    shap_imp_other = pd.DataFrame({
        "feature": feature_cols,
        "shap_other": mean_abs_other
    }).sort_values("shap_other", ascending=False)
    top15_other = shap_imp_other["feature"].head(15).tolist()
    selected = sorted(set(selected + top15_other))
    shap_imp_other.to_csv(os.path.join(output_dir, "shap_other_disease.csv"), index=False)
else:
    print("No Group 2/3 samples.")

print(f"Final selected features: {len(selected)}")

# Save feature lists
gain_imp.to_csv(os.path.join(output_dir, "gain_importance.csv"), index=False)
shap_imp.to_csv(os.path.join(output_dir, "shap_importance.csv"), index=False)
pd.Series(selected).to_csv(os.path.join(output_dir, "selected_features.csv"), index=False, header=["feature"])

# ================== 4. Hyperparameter tuning (5-fold CV, AUC) ==================
print("Hyperparameter tuning (RandomizedSearchCV, 5-fold)...")
X_sel = X[selected]

param_dist = {
    'n_estimators': [50, 100, 200, 300, 400, 500],
    'max_depth': list(range(3, 16)),
    'min_samples_split': list(range(2, 21)),
    'min_samples_leaf': list(range(1, 11)),
    'max_features': np.arange(0.1, 1.0, 0.1)
}

base_rf = RandomForestClassifier(
    random_state=42, class_weight='balanced_subsample',
    bootstrap=True, n_jobs=10
)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
search = RandomizedSearchCV(
    base_rf, param_distributions=param_dist, n_iter=200,
    scoring='roc_auc', cv=cv, verbose=1, random_state=42, n_jobs=10
)
search.fit(X_sel, y)

best_params = search.best_params_
print("Best parameters:", best_params)
print("Best CV AUC:", search.best_score_)

# ================== 5. Final model on all ZMD+GZ ==================
final_model = RandomForestClassifier(
    **best_params,
    random_state=42,
    class_weight='balanced_subsample',
    bootstrap=True,
    n_jobs=10
)
final_model.fit(X_sel, y)

# Save model and feature list
joblib.dump(final_model, os.path.join(output_dir, "final_model.pkl"))
pd.Series(selected).to_csv(os.path.join(output_dir, "model_features.csv"), index=False, header=["feature"])
print(f"Model and features saved to {output_dir}")