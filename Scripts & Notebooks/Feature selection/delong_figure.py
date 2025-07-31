import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_nb_f(mydf):
    p_lst = mydf['Delong2'].tolist()
    i = 0
    while (i + 1 < len(p_lst)) and ((p_lst[i] < 0.05) or (p_lst[i + 1] < 0.05)):
        i += 1
    return i

def main():
    result_path = './result/'  # 默认目录，可根据需要修改

    imp_df = pd.read_csv(os.path.join(result_path, 'Top_90_InfoGain_features.csv'), usecols=['Analytes', 'Ensemble_cv'])
    imp_df.rename(columns={'Ensemble_cv': 'sRNA_imp'}, inplace=True)

    auc_df = pd.read_csv(os.path.join(result_path, 'Delong_Selection_Results.csv'))
    mydf = pd.merge(auc_df, imp_df, how='left', on='Analytes')

    mydf['AUC_lower'] = mydf['AUC_mean'] - mydf['AUC_std']
    mydf['AUC_upper'] = mydf['AUC_mean'] + mydf['AUC_std']
    mydf['AUC_upper'] = mydf['AUC_upper'].clip(upper=1)
    mydf['rna_idx'] = np.arange(1, len(mydf) + 1)

    nb_f = get_nb_f(mydf)

    fig, ax = plt.subplots(figsize=(18, 6.5))
    palette = sns.color_palette("Blues", n_colors=len(mydf))
    palette.reverse()

    sns.barplot(ax=ax, x="Analytes", y="sRNA_imp", palette=palette, data=mydf.sort_values(by="sRNA_imp", ascending=False))

    y_imp_up_lim = round(mydf['sRNA_imp'].max() + 0.01, 2)
    ax.set_ylim([0, y_imp_up_lim])
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xticklabels(mydf['Analytes'], rotation=45, fontsize=10, horizontalalignment='right')

    my_col = ['r'] * nb_f + ['k'] * (len(mydf) - nb_f)
    for ticklabel, tickcolor in zip(ax.get_xticklabels(), my_col):
        ticklabel.set_color(tickcolor)

    ax.set_ylabel('sRNA importance', weight='bold', fontsize=18)
    ax.set_xlabel('')
    ax.grid(which='minor', alpha=0.2, linestyle=':')
    ax.grid(which='major', alpha=0.5, linestyle='--')
    ax.set_axisbelow(True)

    ax2 = ax.twinx()
    ax2.plot(np.arange(nb_f + 1), mydf['AUC_mean'][:nb_f + 1], 'red', alpha=0.8, marker='o')
    ax2.plot(np.arange(nb_f + 1, len(mydf)), mydf['AUC_mean'][nb_f + 1:], 'black', alpha=0.8, marker='o')
    ax2.plot([nb_f, nb_f + 1], mydf['AUC_mean'][nb_f:nb_f + 2], 'black', alpha=0.8, marker='o')

    ax2.fill_between(mydf['rna_idx'] - 1, mydf['AUC_lower'], mydf['AUC_upper'], color='tomato', alpha=0.2)
    ax2.set_ylabel('Cumulative AUC', weight='bold', fontsize=18)
    ax2.tick_params(axis='y', labelsize=14)

    y_auc_up_lim = round(mydf['AUC_upper'].max() + 0.01, 2)
    y_auc_low_lim = round(mydf['AUC_lower'].min() - 0.01, 2)
    ax2.set_ylim([y_auc_low_lim, y_auc_up_lim])

    fig.tight_layout()
    plt.xlim([-0.6, len(mydf) - 0.2])

    out_file = os.path.join(result_path, 'Delong_Plot.svg')
    plt.savefig(out_file, dpi=300, format='svg')
    plt.show()

if __name__ == '__main__':
    main()