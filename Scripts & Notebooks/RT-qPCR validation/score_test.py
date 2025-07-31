import os
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, accuracy_score, confusion_matrix,
    precision_score, f1_score, matthews_corrcoef, roc_curve
)
from sklearn.utils import resample


def sensitivity_score(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[1, 1] / (cm[1, 1] + cm[1, 0])

def specificity_score(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[0, 0] / (cm[0, 0] + cm[0, 1])

def youden_index_score(y_true, y_pred):
    sensitivity = sensitivity_score(y_true, y_pred)
    specificity = specificity_score(y_true, y_pred)
    return sensitivity + specificity - 1

result_path = './rsrna1/'
model_dir = os.path.join(result_path, 'model')

model_files = [f for f in os.listdir(model_dir) if f.endswith('.pkl')]
best_models = {
    os.path.splitext(f)[0].replace('_', ' '): joblib.load(os.path.join(model_dir, f))
    for f in model_files
}

metrics_results = []
roc_data = {}
confusion_matrices = {}

bootstrap_iterations = 1000
confidence_level = 95

for model_name, model in best_models.items():
    print(f"\n=== Evaluating {model_name} ===")
    y_pred = model.predict(X_val)
    y_proba = model.predict_proba(X_val)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X_val)

    metrics = {
        'Model': model_name,
        'AUC': roc_auc_score(y_val, y_proba),
        'Accuracy': accuracy_score(y_val, y_pred),
        'Sensitivity': sensitivity_score(y_val, y_pred),
        'Specificity': specificity_score(y_val, y_pred),
        'Youden Index': youden_index_score(y_val, y_pred),
        'PPV': precision_score(y_val, y_pred),
        'NPV': precision_score(y_val, y_pred, pos_label=0),
        'F1 Score': f1_score(y_val, y_pred),
        'MCC': matthews_corrcoef(y_val, y_pred)
    }

    fpr, tpr, _ = roc_curve(y_val, y_proba)
    roc_data[model_name] = (fpr, tpr, metrics['AUC'])
    confusion_matrices[model_name] = confusion_matrix(y_val, y_pred)

    bootstrap_metrics = {key: [] for key in metrics if key != 'Model'}

    for _ in range(bootstrap_iterations):
        indices = resample(np.arange(len(y_val)), replace=True)
        y_val_boot = y_val.iloc[indices]
        y_pred_boot = y_pred[indices]
        y_proba_boot = y_proba[indices]

        try:
            bootstrap_metrics['AUC'].append(roc_auc_score(y_val_boot, y_proba_boot))
        except:
            bootstrap_metrics['AUC'].append(np.nan)
        bootstrap_metrics['Accuracy'].append(accuracy_score(y_val_boot, y_pred_boot))
        bootstrap_metrics['Sensitivity'].append(sensitivity_score(y_val_boot, y_pred_boot))
        bootstrap_metrics['Specificity'].append(specificity_score(y_val_boot, y_pred_boot))
        bootstrap_metrics['Youden Index'].append(youden_index_score(y_val_boot, y_pred_boot))
        bootstrap_metrics['PPV'].append(precision_score(y_val_boot, y_pred_boot))
        bootstrap_metrics['NPV'].append(precision_score(y_val_boot, y_pred_boot, pos_label=0))
        bootstrap_metrics['F1 Score'].append(f1_score(y_val_boot, y_pred_boot))
        bootstrap_metrics['MCC'].append(matthews_corrcoef(y_val_boot, y_pred_boot))

    for metric in bootstrap_metrics:
        lower = np.nanpercentile(bootstrap_metrics[metric], (100 - confidence_level)/2)
        upper = np.nanpercentile(bootstrap_metrics[metric], 100 - (100 - confidence_level)/2)
        metrics[f'{metric}_lower'] = lower
        metrics[f'{metric}_upper'] = upper

    metrics_results.append(metrics)

metrics_df = pd.DataFrame(metrics_results)
metrics_df.to_csv(os.path.join(result_path, 'validation_metrics.csv'), index=False)