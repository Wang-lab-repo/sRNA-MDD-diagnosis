"""
evaluate.py
===========
External validation of the locked classifier.

Pipeline:
  1. Load trained model and feature list
  2. Fit Platt scaling (logistic calibration) on training set probabilities
  3. Determine optimal Youden threshold on training set
  4. For each test cohort and center:
     - Compute calibrated probabilities via Platt scaling
     - Apply locked threshold for binary classification
     - Report Accuracy, Sensitivity, Specificity, PPV, NPV, F1, AUC, AUPRC
     - Bootstrap 95% CIs (1000 iterations)

Input:
  ./result/train_model/final_model.pkl
  ./result/train_model/model_features.csv
  ./data/ZMD_GZ.csv
  ./data/NJBold_LYG.csv

Output:
  ./result/test_evaluation/test_performance.csv
"""
import pandas as pd
import numpy as np
import re
import os
import joblib
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score,
    precision_score, recall_score, f1_score, confusion_matrix
)
from sklearn.linear_model import LogisticRegression

def clean_column_name(name):
    return re.sub(r'[^\w]', '_', str(name))

# ================== Paths ==================
model_dir = "./result/train_model"
model_path = os.path.join(model_dir, "final_model.pkl")
feature_path = os.path.join(model_dir, "model_features.csv")
train_data_path = "./data/ZMD_GZ.csv"
test_data_paths = [
    "./data/NJBold_LYG.csv",
    # "./data/NJBnew_Tower.csv"
]
output_dir = "./result/test_evaluation"
os.makedirs(output_dir, exist_ok=True)

# ================== 1. Load model and features ==================
model = joblib.load(model_path)
selected_features = pd.read_csv(feature_path)['feature'].tolist()
print(f"Model loaded with {len(selected_features)} features")

# ================== 2. Prepare training set for calibration ==================
train_df = pd.read_csv(train_data_path, index_col=0)
train_df.rename(columns=clean_column_name, inplace=True)
y_train = (train_df['Group'] == 1).astype(int).values
X_train = train_df[selected_features]

# ================== 3. Platt scaling on training set ==================
p_train_raw = model.predict_proba(X_train)[:, 1]
eps = 1e-3
p_train_clipped = np.clip(p_train_raw, eps, 1 - eps)
logit_train = np.log(p_train_clipped / (1 - p_train_clipped))
platt_cal = LogisticRegression(penalty='l2', fit_intercept=True, max_iter=2000)
platt_cal.fit(logit_train.reshape(-1, 1), y_train)
p_train_cal = platt_cal.predict_proba(logit_train.reshape(-1, 1))[:, 1]

# ================== 4. Youden threshold on training set ==================
def find_optimal_threshold(y_true, y_prob, thresholds=np.linspace(0.01, 0.99, 99)):
    best_thr, best_youden = 0.5, -1.0
    for thr in thresholds:
        y_pred = (y_prob >= thr).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0,1]).ravel()
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        youden = sens + spec - 1
        if youden > best_youden:
            best_youden = youden
            best_thr = thr
    return best_thr, best_youden

best_threshold, best_youden = find_optimal_threshold(y_train, p_train_cal)
print(f"Training set Youden threshold: {best_threshold:.3f} (Youden = {best_youden:.3f})")

# ================== 5. Evaluation functions ==================
def compute_metrics(y_true, y_pred, y_prob):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0,1]).ravel()
    return {
        'Accuracy': accuracy_score(y_true, y_pred),
        'Sensitivity': tp / (tp + fn) if (tp + fn) > 0 else 0.0,
        'Specificity': tn / (tn + fp) if (tn + fp) > 0 else 0.0,
        'PPV': precision_score(y_true, y_pred, zero_division=0),
        'NPV': tn / (tn + fn) if (tn + fn) > 0 else 0.0,
        'F1': f1_score(y_true, y_pred, zero_division=0),
        'AUC': roc_auc_score(y_true, y_prob),
        'AUPRC': average_precision_score(y_true, y_prob)
    }

def bootstrap_ci(y_true, y_prob, threshold, n_boot=1000, alpha=0.05, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y_true)
    boot = {k: [] for k in ['Accuracy','Sensitivity','Specificity','PPV','NPV','F1','AUC','AUPRC']}
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        yt, yp = y_true[idx], y_prob[idx]
        yc = (yp >= threshold).astype(int)
        if len(np.unique(yt)) >= 2:
            boot['AUC'].append(roc_auc_score(yt, yp))
            boot['AUPRC'].append(average_precision_score(yt, yp))
        else:
            boot['AUC'].append(np.nan)
            boot['AUPRC'].append(np.nan)
        tn, fp, fn, tp = confusion_matrix(yt, yc, labels=[0,1]).ravel()
        boot['Accuracy'].append(accuracy_score(yt, yc))
        boot['Sensitivity'].append(tp/(tp+fn) if (tp+fn)>0 else 0.0)
        boot['Specificity'].append(tn/(tn+fp) if (tn+fp)>0 else 0.0)
        boot['PPV'].append(tp/(tp+fp) if (tp+fp)>0 else 0.0)
        boot['NPV'].append(tn/(tn+fn) if (tn+fn)>0 else 0.0)
        boot['F1'].append(f1_score(yt, yc, zero_division=0))
    ci = {}
    for metric, vals in boot.items():
        arr = np.array([v for v in vals if not np.isnan(v)])
        if len(arr) == 0:
            mean = low = high = np.nan
        else:
            mean = np.mean(arr)
            low = np.percentile(arr, 100*alpha/2)
            high = np.percentile(arr, 100*(1-alpha/2))
        ci[metric] = (low, high)
    return ci

# ================== 6. Evaluate per center ==================
results = []
for path in test_data_paths:
    test_df = pd.read_csv(path, index_col=0)
    test_df.rename(columns=clean_column_name, inplace=True)
    y_all = (test_df['Group'] == 1).astype(int).values
    X_all = test_df[selected_features]
    hospitals = test_df['Hospital'].unique()

    for center in hospitals:
        mask = test_df['Hospital'] == center
        y_test = y_all[mask]
        X_test = X_all[mask]

        # Raw -> Platt calibrated probability
        p_raw = model.predict_proba(X_test)[:, 1]
        p_clip = np.clip(p_raw, eps, 1 - eps)
        logit_raw = np.log(p_clip / (1 - p_clip))
        y_prob = platt_cal.predict_proba(logit_raw.reshape(-1, 1))[:, 1]

        y_pred = (y_prob >= best_threshold).astype(int)

        base = compute_metrics(y_test, y_pred, y_prob)
        ci = bootstrap_ci(y_test, y_prob, best_threshold)

        row = {'Center': center, 'Threshold': best_threshold,
               'N': len(y_test), 'N_pos': int(y_test.sum())}
        for m in ['Accuracy','Sensitivity','Specificity','PPV','NPV','F1','AUC','AUPRC']:
            low, high = ci[m]
            row[m] = f"{base[m]:.3f} ({low:.3f}, {high:.3f})"
        results.append(row)

# Save and print
res_df = pd.DataFrame(results)
res_csv = os.path.join(output_dir, "test_performance.csv")
res_df.to_csv(res_csv, index=False)
print(res_df.to_string(index=False))
print(f"\nResults saved to {res_csv}")