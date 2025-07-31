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



mirna = pd.read_csv(r'mirna.csv', encoding='GBK')
train_mirna = mirna[mirna.iloc[:, 0].str.startswith('train')]
test_mirna = mirna[mirna.iloc[:, 0].str.startswith('test')]

train_mirna_filtered = train_mirna[train_mirna['group'].isin([0, 1])]
mirna_group = train_mirna_filtered['group']
mirna_feature = train_mirna_filtered.drop(columns=['miRNA', 'group'])

param_grids = {
    'Logistic Regression': {
        'C': [0.01, 0.1, 1, 10, 100],
        'solver': ['saga'],
        'penalty': ['l2'],
        'max_iter': [10000]
    },
    'SVM': {
        'C': [0.1, 1, 10, 100],
        'kernel': ['rbf'],
        'gamma': ['scale', 'auto']
    },
    'Random Forest': {
        'n_estimators': [100, 200, 300],
        'max_depth': [None, 10, 20, 30],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4]
    },
    'MLP': {
        'hidden_layer_sizes': [(10,), (50,), (10, 50), (50, 50)],
        'activation': ['relu', 'tanh'],
        'solver': ['adam'],
        'alpha': [0.0001, 0.001, 0.01],
        'learning_rate': ['constant', 'adaptive']
    },
    'GBDT': {
        'n_estimators': [100, 200, 300],
        'learning_rate': [0.01, 0.1, 0.2],
        'max_depth': [3, 5, 7],
        'subsample': [0.7, 0.8, 1.0]
    },
    'XGBoost': {
        'n_estimators': [100, 200, 300],
        'learning_rate': [0.01, 0.1, 0.2],
        'max_depth': [3, 5, 7],
        'subsample': [0.7, 0.8, 1.0]
    },
    'LightGBM': {
        'n_estimators': [100, 200, 300],
        'learning_rate': [0.01, 0.1, 0.2],
        'num_leaves': [31, 50, 100],
        'boosting_type': ['gbdt'],
        'subsample': [0.7, 0.8, 1.0]
    },
    'CatBoost': {
        'iterations': [100, 200, 300],
        'learning_rate': [0.01, 0.1, 0.2],
        'depth': [3, 5, 7]
    },
    'AdaBoost': {
        'n_estimators': [50, 100, 200],
        'learning_rate': [0.01, 0.1, 1.0]
    }
}

models = {
    'Logistic Regression': LogisticRegression(),
    'SVM': SVC(probability=True),
    'Random Forest': RandomForestClassifier(),
    'MLP': MLPClassifier(),
    'GBDT': GradientBoostingClassifier(),
    'XGBoost': xgb.XGBClassifier(),
    'LightGBM': lgb.LGBMClassifier(),
    'CatBoost': CatBoostClassifier(verbose=0),
    'AdaBoost': AdaBoostClassifier()
}
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
best_params = {}
best_models = {}

for model_name, model in models.items():
    if model_name in param_grids:
        param_grid = param_grids[model_name]
        grid_search = GridSearchCV(model, param_grid, cv=cv, scoring='roc_auc', n_jobs=-1)
        grid_search.fit(mirna_feature, mirna_group)
        best_params[model_name] = grid_search.best_params_
        best_models[model_name] = grid_search.best_estimator_
        print(f"Best parameters for {model_name}: {grid_search.best_params_}")
    else:
        model.fit(mirna_feature, mirna_group)
        best_models[model_name] = model
        print(f"{model_name} has no hyperparameters to tune.")
        

best_params_df = pd.DataFrame(list(best_params.items()), columns=['Model', 'Best Parameters'])
best_params_df['Best Parameters'] = best_params_df['Best Parameters'].apply(lambda x: str(x))  # 转换为字符串
best_params_df.to_csv('./best_params.csv', index=False, encoding='utf-8-sig')

print("Best parameters saved to 'best_params.csv'")



from sklearn.utils import resample
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, f1_score, matthews_corrcoef, roc_curve, confusion_matrix
from sklearn.model_selection import StratifiedKFold
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

results = []
roc_curves = {}

for name, model in best_models.items():
    auc_scores = []
    accuracy_scores = []
    sensitivity_scores = []
    specificity_scores = []
    youden_index_scores = []
    ppv_scores = []
    npv_scores = []
    f1_scores = []
    mcc_scores = []

    y_true_all = []
    y_pred_proba_all = []

    for train_index, test_index in cv.split(mirna_feature, mirna_group):
        X_train_fold, X_test_fold = mirna_feature.iloc[train_index], mirna_feature.iloc[test_index]
        y_train_fold, y_test_fold = mirna_group.iloc[train_index], mirna_group.iloc[test_index]

        model.fit(X_train_fold, y_train_fold)
        y_pred_fold = model.predict(X_test_fold)
        y_pred_proba_fold = model.predict_proba(X_test_fold)[:, 1]

        y_true_all.extend(y_test_fold)
        y_pred_proba_all.extend(y_pred_proba_fold)

        auc_scores.append(roc_auc_score(y_test_fold, y_pred_proba_fold))
        accuracy_scores.append(accuracy_score(y_test_fold, y_pred_fold))
        sensitivity_scores.append(sensitivity_score(y_test_fold, y_pred_fold))
        specificity_scores.append(specificity_score(y_test_fold, y_pred_fold))
        youden_index_scores.append(youden_index_score(y_test_fold, y_pred_fold))
        ppv_scores.append(precision_score(y_test_fold, y_pred_fold))
        npv_scores.append(precision_score(y_test_fold, y_pred_fold, pos_label=0))
        f1_scores.append(f1_score(y_test_fold, y_pred_fold))
        mcc_scores.append(matthews_corrcoef(y_test_fold, y_pred_fold))

        results.append({
            'Model': name,
            'Fold': len(auc_scores),
            'AUC': auc_scores[-1],
            'Accuracy': accuracy_scores[-1],
            'Sensitivity': sensitivity_scores[-1],
            'Specificity': specificity_scores[-1],
            'Youden Index': youden_index_scores[-1],
            'PPV': ppv_scores[-1],
            'NPV': npv_scores[-1],
            'F1 Score': f1_scores[-1],
            'MCC': mcc_scores[-1]
        })

    fpr, tpr, _ = roc_curve(y_true_all, y_pred_proba_all)
    roc_curves[name] = (fpr, tpr, roc_auc_score(y_true_all, y_pred_proba_all))

bootstrap_iterations = 1000
confidence_level = 95
bootstrap_results = {name: [] for name in best_models.keys()}
conf_intervals = {}

for name, model in best_models.items():
    y_true_all = []
    y_pred_proba_all = []

    for train_index, test_index in cv.split(mirna_feature, mirna_group):
        X_train_fold, X_test_fold = mirna_feature.iloc[train_index], mirna_feature.iloc[test_index]
        y_train_fold, y_test_fold = mirna_group.iloc[train_index], mirna_group.iloc[test_index]

        model.fit(X_train_fold, y_train_fold)
        y_pred_proba_fold = model.predict_proba(X_test_fold)[:, 1]

        y_true_all.extend(y_test_fold)
        y_pred_proba_all.extend(y_pred_proba_fold)

    y_true_all = np.array(y_true_all)
    y_pred_proba_all = np.array(y_pred_proba_all)

    for _ in range(bootstrap_iterations):
        y_true_bootstrap, y_pred_proba_bootstrap = resample(y_true_all, y_pred_proba_all)
        bootstrap_results[name].append(roc_auc_score(y_true_bootstrap, y_pred_proba_bootstrap))

    lower_bound = np.percentile(bootstrap_results[name], (100 - confidence_level) / 2)
    upper_bound = np.percentile(bootstrap_results[name], 100 - (100 - confidence_level) / 2)

    conf_intervals[name] = (lower_bound, upper_bound)

    print(f"{name} AUC 95% CI: {lower_bound:.4f} - {upper_bound:.4f}")

metrics = ['Accuracy', 'Sensitivity', 'Specificity', 'Youden Index', 'PPV', 'NPV', 'F1 Score', 'MCC']
bootstrap_metrics_results = {metric: {name: [] for name in best_models.keys()} for metric in metrics}
metrics_conf_intervals = {metric: {} for metric in metrics}

for name, model in best_models.items():
    y_true_all = []
    y_pred_all = []

    for train_index, test_index in cv.split(mirna_feature, mirna_group):
        X_train_fold, X_test_fold = mirna_feature.iloc[train_index], mirna_feature.iloc[test_index]
        y_train_fold, y_test_fold = mirna_group.iloc[train_index], mirna_group.iloc[test_index]

        model.fit(X_train_fold, y_train_fold)
        y_pred_fold = model.predict(X_test_fold)

        y_true_all.extend(y_test_fold)
        y_pred_all.extend(y_pred_fold)

    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)

    for _ in range(bootstrap_iterations):
        y_true_bootstrap, y_pred_bootstrap = resample(y_true_all, y_pred_all)
        
        bootstrap_metrics_results['Accuracy'][name].append(accuracy_score(y_true_bootstrap, y_pred_bootstrap))
        bootstrap_metrics_results['Sensitivity'][name].append(sensitivity_score(y_true_bootstrap, y_pred_bootstrap))
        bootstrap_metrics_results['Specificity'][name].append(specificity_score(y_true_bootstrap, y_pred_bootstrap))
        bootstrap_metrics_results['Youden Index'][name].append(youden_index_score(y_true_bootstrap, y_pred_bootstrap))
        bootstrap_metrics_results['PPV'][name].append(precision_score(y_true_bootstrap, y_pred_bootstrap))
        bootstrap_metrics_results['NPV'][name].append(precision_score(y_true_bootstrap, y_pred_bootstrap, pos_label=0))
        bootstrap_metrics_results['F1 Score'][name].append(f1_score(y_true_bootstrap, y_pred_bootstrap))
        bootstrap_metrics_results['MCC'][name].append(matthews_corrcoef(y_true_bootstrap, y_pred_bootstrap))

    for metric in metrics:
        lower_bound = np.percentile(bootstrap_metrics_results[metric][name], (100 - confidence_level) / 2)
        upper_bound = np.percentile(bootstrap_metrics_results[metric][name], 100 - (100 - confidence_level) / 2)
        metrics_conf_intervals[metric][name] = (lower_bound, upper_bound)
        print(f"{name} {metric} 95% CI: {lower_bound:.4f} - {upper_bound:.4f}")

df_results = pd.DataFrame(results)
df_results.to_csv('./cv_allrna.csv', index=False)

print("Cross-validation results saved to 'cv_results.csv'")

metrics_conf_intervals_df = []
for metric in metrics:
    for name, (lower_bound, upper_bound) in metrics_conf_intervals[metric].items():
        metrics_conf_intervals_df.append({
            'Model': name,
            'Metric': metric,
            'Lower Bound': lower_bound,
            'Upper Bound': upper_bound
        })

df_metrics_conf_intervals = pd.DataFrame(metrics_conf_intervals_df)
df_metrics_conf_intervals.to_csv('./metrics_conf_intervals.csv', index=False)
print("Metrics confidence intervals saved to 'metrics_conf_intervals.csv'")