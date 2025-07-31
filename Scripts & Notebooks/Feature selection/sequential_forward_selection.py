import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import scipy.stats

result_path = 'result/'

top_features_df = pd.read_csv(result_path + 'Top_90_InfoGain_features.csv')
top_features = top_features_df['Analytes'].tolist()
feature_90 = train_df[top_features]


params = {'max_depth': 30, 'min_samples_leaf': 1, 'min_samples_split': 2, 'n_estimators': 100}
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)


def compute_midrank(x):
    """Computes midranks."""
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1)
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T + 1
    return T2


def fastDeLong(predictions_sorted_transposed, label_1_count):
    """
    The fast version of DeLong's method for computing the covariance of
    unadjusted AUC.
    Returns:
       (AUC value, DeLong covariance)
    Reference:
     @article{Guo Y, Chen SD, You J, Huang SY, Chen YL, Zhang Y, Wang LB, He XY, Deng YT, Zhang YR, Huang YY, Dong Q, Feng JF, Cheng W, Yu JT. 
     Multiplex cerebrospinal fluid proteomics identifies biomarkers for diagnosis and prediction of Alzheimer's disease. 
     Nat Hum Behav. 2024 Oct;8(10):2047-2066. doi: 10.1038/s41562-024-01924-6IF: 15.9 Q1 . Epub 2024 Jul 10. PMID: 38987357.
     }
    """

    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty([k, m], dtype=float)
    ty = np.empty([k, n], dtype=float)
    tz = np.empty([k, m + n], dtype=float)
    for r in range(k):
        tx[r, :] = compute_midrank(positive_examples[r, :])
        ty[r, :] = compute_midrank(negative_examples[r, :])
        tz[r, :] = compute_midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def calc_pvalue(aucs, covar):
    """Computes log(10) of p-values.
    Args:
       aucs: 1D array of AUCs
       covar: AUC DeLong covariances
    Returns:
       log10(pvalue)
    """
    l = np.array([[1, -1]])
    z = np.abs(np.diff(aucs)) / np.sqrt(np.dot(np.dot(l, covar), l.T)) 
    return np.log10(2) + scipy.stats.norm.logsf(z, loc=0, scale=1) / np.log(10)


def compute_ground_truth_statistics(ground_truth):
    assert np.array_equal(np.unique(ground_truth), [0, 1])
    order = (-ground_truth).argsort()
    label_1_count = int(ground_truth.sum())
    return order, label_1_count


def delong_roc_variance(ground_truth, predictions):
    """
    Computes ROC AUC variance for a single set of predictions
    Args:
       ground_truth: np.array of 0 and 1
       predictions: np.array of floats of the probability of being class 1
    """
    order, label_1_count = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = predictions[np.newaxis, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    assert len(aucs) == 1
    return aucs[0], delongcov


def delong_roc_test(ground_truth, predictions_one, predictions_two):
    """
    Computes log(p-value) for hypothesis that two ROC AUCs are different
    Args:
       ground_truth: np.array of 0 and 1
       predictions_one: predictions of the first model,
          np.array of floats of the probability of being class 1
       predictions_two: predictions of the second model,
          np.array of floats of the probability of being class 1
    """
    order, label_1_count = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = np.vstack((predictions_one, predictions_two))[:, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    return calc_pvalue(aucs, delongcov)


y_pred_prev1 = np.zeros(len(df_group))
y_pred_prev2 = np.zeros(len(df_group))
y_pred_prev3 = np.zeros(len(df_group))
selected_features = []
results = []

for f in top_features:
    selected_features.append(f)
    X_subset = feature_90[selected_features]

    y_pred_all = []
    y_true_all = []
    auc_fold = []

    for train_idx, test_idx in cv.split(X_subset, df_group):
        X_train, X_test = X_subset.iloc[train_idx], X_subset.iloc[test_idx]
        y_train, y_test = df_group.iloc[train_idx], df_group.iloc[test_idx]

        model = RandomForestClassifier(**params)
        model.fit(X_train, y_train)
        y_pred_prob = model.predict_proba(X_test)[:, 1]

        auc_fold.append(roc_auc_score(y_test, y_pred_prob))
        y_pred_all += y_pred_prob.tolist()
        y_true_all += y_test.tolist()

    auc_all = roc_auc_score(y_true_all, y_pred_all)
    log10_p1 = delong_roc_test(np.array(y_true_all), np.array(y_pred_prev1), np.array(y_pred_all))
    log10_p2 = delong_roc_test(np.array(y_true_all), np.array(y_pred_prev2), np.array(y_pred_all))
    log10_p3 = delong_roc_test(np.array(y_true_all), np.array(y_pred_prev3), np.array(y_pred_all))

    y_pred_prev3 = y_pred_prev2
    y_pred_prev2 = y_pred_prev1
    y_pred_prev1 = y_pred_all

    out = np.array([
        np.mean(auc_fold),
        np.std(auc_fold),
        10**log10_p1[0][0],
        10**log10_p2[0][0],
        10**log10_p3[0][0],
        auc_all
    ])
    results.append(out)

    print(f"Feature: {f}, AUC: {out}")

result_df = pd.DataFrame(results, columns=['AUC_mean', 'AUC_std', 'Delong1', 'Delong2', 'Delong3', 'AUC_all'])
result_df = pd.concat([pd.DataFrame({'Analytes': selected_features}), result_df], axis=1)
result_df.to_csv(result_path + 'Delong_Selection_Results.csv', index=False)

imp_df = pd.read_csv(result_path + 'Top_90_InfoGain_features.csv', usecols=['Analytes', 'Ensemble_cv'])
imp_df.rename(columns={'Ensemble_cv': 'sRNA_imp'}, inplace=True)
merged_df = pd.merge(result_df, imp_df, on='Analytes', how='left')

merged_df.to_csv(result_path + 'Merged_AUC_Importance.csv', index=False)

print("Finished feature selection and merged with importance.")