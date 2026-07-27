"""
plot_figures.py
===============
Generate ROC, Decision Curve Analysis (DCA), and Calibration plots.

For each test cohort and center, produces:
  - ROC curve with AUC
  - DCA curve (net benefit vs. threshold)
  - Calibration curve with Brier score, slope, and intercept

Requires the trained model and Platt scaling to be re-fitted on training data.

Input:
  ./result/train_model/final_model.pkl
  ./result/train_model/model_features.csv
  ./data/ZMD_GZ.csv
  ./data/NJBold_LYG.csv

Output:
  ./result/figures/ROC_*.pdf
  ./result/figures/DCA_*.pdf
  ./result/figures/Calibration_*.pdf
"""
import pandas as pd
import numpy as np
import re
import os
import joblib
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, brier_score_loss
from sklearn.calibration import calibration_curve

# ================== Global plot settings ==================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
})
DARK_GREEN = "#1D5A78"
RED = "#9D1E1E"

def clean_column_name(name):
    return re.sub(r'[^\w]', '_', str(name))

# ================== Paths ==================
model_path = './result/train_model/final_model.pkl'
feature_path = './result/train_model/model_features.csv'
train_data_path = './data/ZMD_GZ.csv'
test_files = {
    'NJBold_LYG': './data/NJBold_LYG.csv',
    # 'NJBnew_Tower': './data/NJBnew_Tower.csv'
}
output_dir = './result/figures'
os.makedirs(output_dir, exist_ok=True)

# ================== Load model and features ==================
model = joblib.load(model_path)
selected_features = pd.read_csv(feature_path)['feature'].tolist()

# ================== Fit Platt scaling on training set ==================
train_df = pd.read_csv(train_data_path, index_col=0)
train_df.rename(columns=clean_column_name, inplace=True)
y_train = (train_df['Group'] == 1).astype(int).values
X_train = train_df[selected_features]
p_train_raw = model.predict_proba(X_train)[:, 1]
eps = 1e-3
p_train_clip = np.clip(p_train_raw, eps, 1 - eps)
logit_train = np.log(p_train_clip / (1 - p_train_clip))
platt = LogisticRegression(penalty='l2', fit_intercept=True, max_iter=2000)
platt.fit(logit_train.reshape(-1, 1), y_train)

# ================== Functions for net benefit and DCA ==================
def net_benefit(y_true, y_prob, threshold):
    n = len(y_true)
    y_pred = (y_prob >= threshold).astype(int)
    tp = np.sum((y_pred == 1) & (y_true == 1))
    fp = np.sum((y_pred == 1) & (y_true == 0))
    nb = (tp - fp * threshold / (1 - threshold)) / n
    return nb

def dca_curve(y_true, y_prob, thresholds):
    nb_model = np.zeros_like(thresholds)
    nb_all = np.zeros_like(thresholds)
    prev = y_true.mean()
    for i, pt in enumerate(thresholds):
        if pt == 0:
            nb_model[i] = prev
            nb_all[i] = prev
        else:
            nb_model[i] = net_benefit(y_true, y_prob, pt)
            nb_all[i] = prev - (1 - prev) * pt / (1 - pt)
    return nb_model, nb_all

# ================== Process each center in each test file ==================
for cohort_name, data_path in test_files.items():
    test_df = pd.read_csv(data_path, index_col=0)
    test_df.rename(columns=clean_column_name, inplace=True)
    y_all = (test_df['Group'] == 1).astype(int)
    X_all = test_df[selected_features]

    # Split by Hospital (center)
    hospitals = test_df['Hospital'].unique()
    for center in hospitals:
        print(f"Processing {center}...")
        mask = test_df['Hospital'] == center
        y_test = y_all[mask].values
        X_test = X_all[mask]

        # Raw probabilities -> Platt calibrated
        p_raw = model.predict_proba(X_test)[:, 1]
        p_clip = np.clip(p_raw, eps, 1 - eps)
        logit = np.log(p_clip / (1 - p_clip))
        y_prob = platt.predict_proba(logit.reshape(-1, 1))[:, 1]

        # ---- ROC curve ----
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        roc_auc = auc(fpr, tpr)
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.plot(fpr, tpr, color=DARK_GREEN, lw=1.5,
                label=f'{center} (AUC = {roc_auc:.3f})')
        ax.plot([0, 1], [0, 1], color='grey', lw=0.8, linestyle='--', alpha=0.7)
        ax.set_xlim([0.0, 1.02])
        ax.set_ylim([0.0, 1.02])
        ax.set_xlabel('1 - Specificity')
        ax.set_ylabel('Sensitivity')
        ax.set_title(f'ROC - {center}', fontweight='bold')
        ax.legend(loc='lower right', frameon=False)
        ax.grid(axis='both', linestyle='--', linewidth=0.3, alpha=0.4)
        plt.tight_layout()
        roc_path = os.path.join(output_dir, f'ROC_{center}.pdf')
        fig.savefig(roc_path)
        plt.close(fig)

        # ---- DCA ----
        thresholds = np.linspace(0.01, 0.99, 100)
        nb_model, nb_all = dca_curve(y_test, y_prob, thresholds)
        fig, ax = plt.subplots(figsize=(4.5, 4))
        ax.plot(thresholds, nb_model, color=DARK_GREEN, lw=1.5, label='Model')
        ax.plot(thresholds, nb_all, color=RED, linestyle='--', lw=1.0, label='Treat All')
        ax.plot(thresholds, np.zeros_like(thresholds), color='grey', linestyle=':', lw=1.0, label='Treat None')
        prob_min = np.percentile(y_prob, 1)
        prob_max = np.percentile(y_prob, 99)
        ax.set_xlim(max(0.0, prob_min - 0.02), min(1.0, prob_max + 0.02))
        ymax = max(nb_model.max(), nb_all.max()) + 0.05
        ax.set_ylim(-0.1, ymax)
        ax.set_xlabel('Threshold Probability')
        ax.set_ylabel('Net Benefit')
        ax.set_title(f'DCA - {center}', fontweight='bold')
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='grey')
        ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.5)
        plt.tight_layout()
        dca_path = os.path.join(output_dir, f'DCA_{center}.pdf')
        fig.savefig(dca_path)
        plt.close(fig)

        # ---- Calibration curve ----
        n_bins = 5
        prob_true, prob_pred = calibration_curve(y_test, y_prob, n_bins=n_bins, strategy='quantile')
        brier = brier_score_loss(y_test, y_prob)
        logit_test = np.log(np.clip(y_prob, eps, 1 - eps) / (1 - np.clip(y_prob, eps, 1 - eps)))
        lr_cal = LogisticRegression(penalty=None, fit_intercept=True, max_iter=2000)
        lr_cal.fit(logit_test.reshape(-1, 1), y_test)
        slope = lr_cal.coef_[0][0]
        intercept = lr_cal.intercept_[0]

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.plot([0, 1], [0, 1], 'k--', lw=1.0, label='Ideal')
        ax.plot(prob_pred, prob_true, marker='o', color=DARK_GREEN,
                lw=1.5, markersize=6, markerfacecolor='white', label='Calibrated')
        text_str = f'Brier = {brier:.3f}\nSlope = {slope:.3f}\nIntercept = {intercept:.3f}'
        ax.text(0.05, 0.85, text_str, transform=ax.transAxes, fontsize=7,
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='grey', alpha=0.8))
        ax.set_title(f'Calibration - {center}', fontweight='bold')
        ax.set_xlabel('Predicted Probability')
        ax.set_ylabel('Observed Proportion')
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)
        ax.set_aspect('equal')
        ax.legend(loc='lower right', frameon=True, facecolor='white', edgecolor='grey')
        plt.tight_layout()
        cal_path = os.path.join(output_dir, f'Calibration_{center}.pdf')
        fig.savefig(cal_path)
        plt.close(fig)

        print(f"{center}: ROC={roc_auc:.3f}, Brier={brier:.3f}")

print(f"All figures saved to {output_dir}")