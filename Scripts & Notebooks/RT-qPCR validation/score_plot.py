import os
import pandas as pd
import numpy as np
import joblib
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

result_path = './rsrna1/'
model_dir = os.path.join(result_path, 'model')
csv_file = 'rsrna1.csv'
miRNA_list = ['rsrna1']
encoding = 'GBK'
groups = ['MDD', 'Control']
target_group = "MDD"
palette = {
    'MDD': '#FFA500',
    'Control': '#2E8B57'
}


# Select the best model based on NSGA optimization
best_model = best_models[best_model_name]

y_proba = best_model.predict_proba(X_val)[:, 1]

prob_df = pd.DataFrame({
    'True_label': y_val,
    'Predicted_probability': y_proba
})
prob_path = os.path.join(result_path, 'sample_probabilities.csv')
prob_df.to_csv(prob_path, index=False)

all = prob_df.copy()
selected_mirnas = ['Predicted_probability']

for mirna in selected_mirnas:
    value_df = all.dropna(subset=[mirna, 'True_label']).copy()
    value_df = value_df[value_df['True_label'].isin(groups)]

    filtered_data = pd.DataFrame()
    for group in groups:
        group_data = value_df[value_df['True_label'] == group]
        Q1 = group_data[mirna].quantile(0.25)
        Q3 = group_data[mirna].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        filtered_group = group_data[(group_data[mirna] >= lower_bound) & (group_data[mirna] <= upper_bound)]
        filtered_data = pd.concat([filtered_data, filtered_group])

    fig, ax = plt.subplots(figsize=(4.5, 6))
    sns.boxplot(
        x='True_label', y=mirna, data=filtered_data, order=groups,
        palette=palette, ax=ax, linewidth=1.5, width=0.35, showfliers=False
    )
    sns.stripplot(
        x='True_label', y=mirna, data=filtered_data, order=groups,
        palette=palette, edgecolor='w', linewidth=0.4,
        size=5, alpha=0.75, jitter=0.15, ax=ax
    )

    ax.tick_params(axis='both', which='both', direction='out', length=6,
                   width=1.5, colors='black', bottom=True, left=True)
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1line.set_visible(True)
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_visible(True)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_title('Predicted Probability Distribution', fontsize=16)
    ax.set_ylabel('Predicted probability', fontsize=12)
    ax.set_xlabel('')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1)
    ax.grid(True, linestyle='--', alpha=0.7)

    compare_groups = [g for g in groups if g != target_group]
    combinations = [(target_group, other) for other in compare_groups]
    y_min, y_max = ax.get_ylim()
    vertical_step = (y_max - y_min) * 0.08

    for i, (group1, group2) in enumerate(combinations):
        data1 = value_df[value_df['True_label'] == group1][mirna]
        data2 = value_df[value_df['True_label'] == group2][mirna]
        levene_stat, levene_p = stats.levene(data1, data2)
        equal_var = levene_p >= 0.05
        t_stat, p_val = stats.ttest_ind(data1, data2, equal_var=equal_var)

        x1 = groups.index(group1)
        x2 = groups.index(group2)
        y_pos = y_max + (i+1)*vertical_step

        ax.plot([x1, x2], [y_pos, y_pos], lw=1.5, color='black')
        if p_val < 0.001:
            p_text = '***'
        elif p_val < 0.01:
            p_text = '**'
        elif p_val < 0.05:
            p_text = '*'
        else:
            p_text = f'p={p_val:.3f}'

        ax.text((x1+x2)/2, y_pos + vertical_step/4, p_text,
                ha='center', va='bottom', fontsize=12,
                backgroundcolor=(1, 1, 1, 0.5))

    ax.set_ylim(y_min, y_max + len(combinations)*vertical_step + vertical_step)

    plt.tight_layout()
    fig_path = f'mdd_{mirna}.svg'
    plt.savefig(fig_path, dpi=600, bbox_inches='tight')
    plt.close()