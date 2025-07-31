from sklearn.metrics import (
    roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, 
    confusion_matrix, average_precision_score, log_loss, matthews_corrcoef
)
import numpy as np
import pandas as pd
import time

def specificity_score(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    return cm[0, 0] / (cm[0, 0] + cm[0, 1])

def evaluate_model_performance(models_dict, dataset_dict):
    results = []

    for name, models in models_dict.items():
        test_features, test_labels = dataset_dict[name]

        for model_name, model in models.items():
            try:
                start = time.time()
                y_pred_prob = model.predict_proba(test_features)[:, 1]
                pred_time = (time.time() - start) / len(test_features)

                y_pred = model.predict(test_features)

                result = {
                    "Dataset": name,
                    "Model": model_name,
                    "AUC": roc_auc_score(test_labels, y_pred_prob),
                    "Accuracy": accuracy_score(test_labels, y_pred),
                    "Sensitivity (SN)": recall_score(test_labels, y_pred),
                    "Specificity (SP)": specificity_score(test_labels, y_pred),
                    "F1 Score": f1_score(test_labels, y_pred),
                    "AUPRC": average_precision_score(test_labels, y_pred_prob),
                    "Log-Loss": log_loss(test_labels, y_pred_prob),
                    "MCC": matthews_corrcoef(test_labels, y_pred),
                    "InferenceTime": pred_time
                }

                results.append(result)
            except Exception as e:
                print(f"Error evaluating {model_name} on {name}: {e}")

    return pd.DataFrame(results)
