import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# take miRNA for example
data = pd.read_csv('mirna.csv', index_col=0)
hamd_scores = data.iloc[0]
mirna_data = data.iloc[1:]

correlations = {}
p_values = {}

for mirna in mirna_data.index:
    correlation, p_value = spearmanr(mirna_data.loc[mirna], hamd_scores)
    correlations[mirna] = correlation
    p_values[mirna] = p_value

result = pd.DataFrame({'Correlation': correlations, 'P_value': p_values})
result['Adjusted_p_value'] = multipletests(result['P_value'], method='fdr_bh')[1]
result.to_csv('./mirna/spearman_results_corrected.csv')

significant_positive = result[(result['Adjusted_p_value'] < 0.05) & (result['Correlation'] > 0.2)]
significant_negative = result[(result['Adjusted_p_value'] < 0.05) & (result['Correlation'] < -0.2)]
combined_significant = pd.concat([significant_positive, significant_negative], axis=0)
combined_significant.to_csv('./mirna/combined_significant_mirna.csv')
