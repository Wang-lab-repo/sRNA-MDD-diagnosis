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
