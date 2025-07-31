import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, accuracy_score, confusion_matrix, roc_curve, 
    precision_score, recall_score, f1_score, matthews_corrcoef
)
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
import xgboost as xgb
import lightgbm as lgb
from catboost import CatBoostClassifier

mirna=pd.read_csv(r'mirna.csv',encoding='GBK')
train_mirna = mirna[mirna.iloc[:, 0].str.startswith('train')]
test_mirna = mirna[mirna.iloc[:, 0].str.startswith('test')]
test_mirna_filtered = test_mirna[test_mirna['group'].isin([0, 1])]
X_test = test_mirna_filtered.drop(columns=['miRNA','group'])
y_test = test_mirna_filtered['group']

from sklearn.utils import resample
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, recall_score, confusion_matrix, precision_score, f1_score, matthews_corrcoef
import pandas as pd

results = []

plt.figure(figsize=(10, 8))

for name, model in best_models.items():
    if hasattr(model, "decision_function"):
        y_score = model.decision_function(X_test)
    else:
        y_score = model.predict_proba(X_test)[:, 1]

    fpr, tpr, thresholds = roc_curve(y_test, y_score)
    auc = roc_auc_score(y_test, y_score)

    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    sensitivity = recall_score(y_test, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    specificity = tn / (tn + fp)
    ppv = precision_score(y_test, y_pred)
    npv = tn / (tn + fn)
    f1 = f1_score(y_test, y_pred)
    mcc = matthews_corrcoef(y_test, y_pred)
    youden_index = sensitivity + specificity - 1
   
    n_bootstraps = 1000
    auc_scores, accuracy_scores, sensitivity_scores, specificity_scores = [], [], [], []
    ppv_scores, npv_scores, f1_scores, mcc_scores, youden_index_scores = [], [], [], [], []
    
    for _ in range(n_bootstraps):
        y_test_resampled, y_score_resampled = resample(y_test, y_score, random_state=np.random.randint(1, 100))
        auc_resampled = roc_auc_score(y_test_resampled, y_score_resampled)
        auc_scores.append(auc_resampled)

        y_pred_resampled = np.where(y_score_resampled > 0.5, 1, 0)
        accuracy_scores.append(accuracy_score(y_test_resampled, y_pred_resampled))
        sensitivity_scores.append(recall_score(y_test_resampled, y_pred_resampled))
        tn_resampled, fp_resampled, fn_resampled, tp_resampled = confusion_matrix(y_test_resampled, y_pred_resampled).ravel()
        specificity_scores.append(tn_resampled / (tn_resampled + fp_resampled))
        ppv_scores.append(precision_score(y_test_resampled, y_pred_resampled))
        npv_scores.append(tn_resampled / (tn_resampled + fn_resampled))
        f1_scores.append(f1_score(y_test_resampled, y_pred_resampled))
        mcc_scores.append(matthews_corrcoef(y_test_resampled, y_pred_resampled))
        youden_index_scores.append(sensitivity_scores[-1] + specificity_scores[-1] - 1)
    
    def compute_ci(scores):
        sorted_scores = np.array(scores)
        sorted_scores.sort()
        lower = sorted_scores[int(0.025 * len(sorted_scores))]
        upper = sorted_scores[int(0.975 * len(sorted_scores))]
        return lower, upper
    
    auc_ci = compute_ci(auc_scores)
    accuracy_ci = compute_ci(accuracy_scores)
    sensitivity_ci = compute_ci(sensitivity_scores)
    specificity_ci = compute_ci(specificity_scores)
    ppv_ci = compute_ci(ppv_scores)
    npv_ci = compute_ci(npv_scores)
    f1_ci = compute_ci(f1_scores)
    mcc_ci = compute_ci(mcc_scores)
    youden_index_ci = compute_ci(youden_index_scores)
    
    plt.plot(fpr, tpr, label=f'{name} (AUC = {auc:.3f} [{auc_ci[0]:.3f}-{auc_ci[1]:.3f}])')

    results.append({
        'Model': name,
        'AUC': f'{auc:.3f} [{auc_ci[0]:.3f}-{auc_ci[1]:.3f}]',
        'Accuracy': f'{accuracy:.3f} [{accuracy_ci[0]:.3f}-{accuracy_ci[1]:.3f}]',
        'Sensitivity': f'{sensitivity:.3f} [{sensitivity_ci[0]:.3f}-{sensitivity_ci[1]:.3f}]',
        'Specificity': f'{specificity:.3f} [{specificity_ci[0]:.3f}-{specificity_ci[1]:.3f}]',
        'Youden Index': f'{youden_index:.3f} [{youden_index_ci[0]:.3f}-{youden_index_ci[1]:.3f}]',
        'PPV': f'{ppv:.3f} [{ppv_ci[0]:.3f}-{ppv_ci[1]:.3f}]',
        'NPV': f'{npv:.3f} [{npv_ci[0]:.3f}-{npv_ci[1]:.3f}]',
        'F1 Score': f'{f1:.3f} [{f1_ci[0]:.3f}-{f1_ci[1]:.3f}]',
        'MCC': f'{mcc:.3f} [{mcc_ci[0]:.3f}-{mcc_ci[1]:.3f}]'
    })

plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.01])

plt.xlabel('False Positive Rate', fontsize=12)
plt.ylabel('True Positive Rate', fontsize=12)
plt.title('ROC Curve', fontsize=15)
plt.legend(loc="lower right")

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.savefig('ROC.svg', format='svg')

plt.show()

results_df = pd.DataFrame(results)
print(results_df)

results_df.to_csv('results.csv', index=False)
