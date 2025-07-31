import os
import numpy as np
import pandas as pd
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
import shap

dpath = ''
result_path = 'result/'

df = pd.read_csv(dpath + 'data/all.csv', encoding='GBK')
train_df = df[df.iloc[:, 0].str.startswith('train')]
train_df_filtered = train_df[train_df['group'].isin([1, 0, 3])].copy()
train_df_filtered['group'] = train_df_filtered['group'].replace({3: 0})
df_group = train_df_filtered['group']
df_feature = train_df_filtered.drop(columns=['allRNA', 'group'])

params = {'max_depth': 30, 'min_samples_leaf': 1, 'min_samples_split': 2, 'n_estimators': 100}
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

def normal_imp(mydict):
    mysum = sum(mydict.values())
    for key in mydict.keys():
        mydict[key] = mydict[key] / mysum
    return mydict

tg_imp_cv = Counter()
shap_imp_cv = np.zeros(df_feature.shape[1])

for train_idx, test_idx in cv.split(df_feature, df_group):
    X_train, X_test = df_feature.iloc[train_idx], df_feature.iloc[test_idx]
    y_train, y_test = df_group.iloc[train_idx], df_group.iloc[test_idx]

    model = RandomForestClassifier(**params)
    model.fit(X_train, y_train)

    gain_importance = dict(zip(df_feature.columns, model.feature_importances_))
    tg_imp_cv += Counter(normal_imp(gain_importance))

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)
    shap_mean = np.mean(np.abs(shap_values[1]), axis=0)
    shap_imp_cv += shap_mean / np.sum(shap_mean)

shap_imp_df = pd.DataFrame({
    'Analytes': df_feature.columns,
    'ShapValues_cv': shap_imp_cv / cv.get_n_splits()
})
shap_imp_df.sort_values(by='ShapValues_cv', ascending=False, inplace=True)

tg_imp_cv = normal_imp(tg_imp_cv)
tg_imp_df = pd.DataFrame({
    'Analytes': list(tg_imp_cv.keys()),
    'TotalGain_cv': list(tg_imp_cv.values())
})

imp_df = pd.merge(shap_imp_df, tg_imp_df, on='Analytes')
imp_df['Ensemble_cv'] = (imp_df['ShapValues_cv'] + imp_df['TotalGain_cv']) / 2
imp_df.sort_values(by='TotalGain_cv', ascending=False, inplace=True)
imp_df.to_csv(result_path + 'mdd_Importance483.csv', index=False)

def get_top_features_by_prop(df, prop=0.9):
    score = 0
    i = 0
    while score < prop and i < len(df):
        score += df['Ensemble_cv'].iloc[i]
        i += 1
    return i

top_n = get_top_features_by_prop(imp_df, top_prop=0.9)
top_features_df = imp_df.iloc[:top_n]
top_features_df.to_csv(result_path + 'Top_90_InfoGain_features.csv', index=False)
