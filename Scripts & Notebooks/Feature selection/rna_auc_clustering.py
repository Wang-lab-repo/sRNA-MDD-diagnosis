import os
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score


def calculate_auc_for_feature(RNA_feature, df_feature, df_group, cv, params):
    aucs = []
    for train_idx, test_idx in cv.split(df_feature, df_group):
        X_train, X_test = df_feature.iloc[train_idx], df_feature.iloc[test_idx]
        y_train, y_test = df_group.iloc[train_idx], df_group.iloc[test_idx]
        
        model = RandomForestClassifier(**params)
        model.fit(X_train[[RNA_feature]], y_train)
        
        y_pred_prob = model.predict_proba(X_test[[RNA_feature]])[:, 1]
        auc = roc_auc_score(y_test, y_pred_prob)
        aucs.append(auc)
    
    return np.mean(aucs)


def main():
    dpath = './'
    result_path = dpath + 'mdd_ranking/'

    df = pd.read_csv(dpath + 'data/all.csv', encoding='GBK')

    train_df = df[df.iloc[:, 0].str.startswith('train')]
    train_df_filtered = train_df[train_df['group'].isin([1, 0, 3])].copy()
    train_df_filtered['group'] = train_df_filtered['group'].replace({3: 0})

    df_group = train_df_filtered['group']
    df_feature = train_df_filtered.drop(columns=['allRNA', 'group'])

    correlation_matrix = df_feature.corr(method='spearman')
    correlation_matrix.to_csv(result_path + 'mdd_correlation_matrix.csv', index=True)

    distance_matrix = 1 - np.abs(correlation_matrix)
    dist_array = squareform(distance_matrix)
    dist_linkage = hierarchy.linkage(dist_array, method='ward')

    cluster_ids = hierarchy.fcluster(dist_linkage, 0.6, criterion="distance")
    feature_names = df_feature.columns.tolist()
    cluster_data = pd.DataFrame({'RNA_feature': feature_names, 'ClusterID': cluster_ids})
    output_csv = result_path + 'mdd_feature_clusters.csv'
    cluster_data.to_csv(output_csv, index=False)

    params = {'max_depth': 30, 'min_samples_leaf': 1, 'min_samples_split': 2, 'n_estimators': 100}
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    auc_scores = {RNA_feature: calculate_auc_for_feature(RNA_feature, df_feature, df_group, cv, params)
                  for RNA_feature in df_feature.columns}

    auc_df = pd.DataFrame(list(auc_scores.items()), columns=['RNA_feature', 'AUC'])
    output_auc_file = result_path + 'mdd_auc_scores.csv'
    auc_df.to_csv(output_auc_file, index=False)
    print(f"AUC scores have been saved to {output_auc_file}")

    feature_auc = pd.read_csv(output_auc_file, encoding='GBK')
    feature_clus = pd.read_csv(output_csv, encoding='GBK')
    feature_sel = pd.merge(feature_auc, feature_clus, how='left', on="RNA_feature")
    best_features = feature_sel.loc[feature_sel.groupby('ClusterID')['AUC'].idxmax()]
    output_best_features = result_path + 'mdd_best_features.csv'
    best_features.to_csv(output_best_features, index=False)


if __name__ == '__main__':
    main()
