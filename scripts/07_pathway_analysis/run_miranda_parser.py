import os
import csv

def TidyMirandaResult(path, inputfile, outfile):
    infpath = os.path.join(path, inputfile)
    outfpath = os.path.join(path, outfile)
    
    with open(infpath, 'r') as f:
        lines = f.readlines()
        
    with open(outfpath, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['sRNAs', 'Gene', 'Score', 'Energy']) 
        
        for line in lines:
            line = line.strip()
            if line.startswith('>>'):
                parts = line.split('\t')

                srna_name = parts[0][2:] 
                gene_name = parts[1].split('|')[0] 
                score = parts[2]
                energy = parts[3] 
                
                csv_writer.writerow([srna_name, gene_name, score, energy])

def main():
    path = r'./'
    input_file = 'miranda_result.txt'
    output_file = 'parsed_miranda.csv'
    TidyMirandaResult(path, input_file, output_file)

if __name__ == "__main__":
    main()
