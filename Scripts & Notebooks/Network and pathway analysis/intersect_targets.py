
import pandas as pd

hybrid_result_df = pd.read_csv('parsed_miranda.csv')
miranda_result_df = pd.read_csv('parsed_rnahybrid.csv')

miranda_result_df = miranda_result_df[miranda_result_df['Energy'] < -20]

sRNAs_list = hybrid_result_df['sRNAs'].unique()

result = []

for srna in sRNAs_list:
    hybrid_genes = set(hybrid_result_df[hybrid_result_df['sRNAs'] == srna]['Gene'])
    
    miranda_genes = set(miranda_result_df[miranda_result_df['sRNAs'] == srna]['Gene'])

    intersect_genes = hybrid_genes.intersection(miranda_genes)

    for gene in intersect_genes:
        result.append([srna, gene])

result_df = pd.DataFrame(result, columns=['sRNA', 'Gene'])

result_df.to_csv('sRNA_gene_intersection.csv', index=False)