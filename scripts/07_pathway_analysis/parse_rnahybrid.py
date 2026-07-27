import re
import pandas as pd

def parse_rnahybrid_results(result_file):
    results = []
    current_srna = None
    current_target = None
    current_pvalue = None
    
    with open(result_file, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.startswith("miRNA :"):
                current_srna = re.search(r'^miRNA : (\S+)', line).group(1)
            
            elif line.startswith("target:"):
                target_line = re.search(r'^target: (\S+)', line).group(1)
                current_target = target_line.split('|')[0]
            
            elif line.startswith("p-value:"):
                current_pvalue = float(re.search(r'p-value: (\S+)', line).group(1))
                
                if current_pvalue and current_pvalue < 0.05:
                    results.append([current_srna, current_target, current_pvalue])
    
    return results

rsRNA_targets = parse_rnahybrid_results('parsed_rnahybrid.csv')

pd.DataFrame(rsRNA_targets, columns=['sRNAs', 'Gene', 'p-value']).drop_duplicates().to_csv("parsed_rnahybrid.csv", index=False)