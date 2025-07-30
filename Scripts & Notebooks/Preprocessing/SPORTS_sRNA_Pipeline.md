
# SPORTS Small RNA Annotation and Quantification Pipeline

This document provides a complete guide for processing and quantifying rsRNA and tsRNA using SPORTS1.1.

---

## Step 1: Setup

```bash
cd /path/to/project
mkdir -p sports pirna mirna

base_dir=/path/to/project
clean_data=$base_dir/clean
sports_data=$base_dir/sports

find $clean_data -name '*.fq'  > $sports_data/sports_name.txt
find $clean_data -name '*.fq' -exec basename {} \; | sed "s/_trimmed.fq//g" | sort > $sports_data/short_sports_name.list
```

---

## Step 2: Run SPORTS1.1

```bash
export PATH=/path/to/sports1.1/source:$PATH

base_dir=/path/to/project
clean_data=$base_dir/clean
genome_dir=/path/to/reference/Mus_musculus
output_dir=$base_dir/sports

sports.pl -i $output_dir/sports_name.txt -p 8 \
-g $genome_dir/genome/genome \
-m $genome_dir/miRBase/miRBase \
-r $genome_dir/rRNAdb/human_rRNA \
-t $genome_dir/GtRNAdb/hg19-tRNAs \
-w $genome_dir/piRBase/piR_human \
-o $output_dir
```

---

## Step 3: Build rsRNA + tsRNA Count Matrix


```bash
cd /path/to/project
mkdir -p sports_merge
find ./sports -type f -name "*R1*output.txt" > R1_sports_name.list
find ./sports -type f -name "*R1*output.txt" -exec basename {} \; > R1_sports_short_name.list

cat R1_sports_name.list | while read id
do
  file_name=$(basename "$id" | cut -d '_' -f 1)
  awk 'BEGIN {OFS=","} {print $6}' "$id" | sort -u | grep -v "NO_Annotation" >> "./sports_merge/all_anno.csv"
done

sort -u "./sports_merge/all_anno.csv" -o "./sports_merge/all_anno.csv"

cat R1_sports_name.list | while read id
do
  file_name=$(basename "$id" | cut -d '_' -f 1)
  awk 'BEGIN {OFS=","} {print $2,$6,$4}' "$id" | grep -v "NO_Annotation" > "./sports_merge/1_${file_name}.csv"
  sed "s/,/:/" "./sports_merge/1_${file_name}.csv" | sed "s/Reads/${file_name}/g" > "./sports_merge/3_${file_name}.csv"
  rm "./sports_merge/1_${file_name}.csv"
done
```

---

## Step 4: Merge Into Final Matrix

```bash
python merge.py
```

---

## Step 5: Subset rsRNA and tsRNA

```bash
cd /path/to/project/sports_R1merge

mkdir -p rsrna tsrna

# rsRNA
input_file=./merged_result.csv
output_file=./rsrna/rsRNA_result.csv
head -n 1 "$input_file" > "$output_file"
awk '/-rRNA/' "$input_file" >> "$output_file"

# tsRNA
output_file=./tsrna/tsRNA_result.csv
head -n 1 "$input_file" > "$output_file"
awk '/tRNA-/' "$input_file" >> "$output_file"
```

---
## Reference Databases
The following reference databases were used for small RNA annotation:  

Genomic tRNA Database (GtRNAdb v21, hg19):  
tRNA gene annotations were retrieved from the Genomic tRNA Database (GtRNAdb), version 21 for the human genome build hg19.  

rRNA Reference Sequences:  
Ribosomal RNA (rRNA) reference sequences were compiled from the NCBI Nucleotide and NCBI Gene databases.  

---

## Notes

- Ensure all reference paths are updated according to your cluster setup.  
- SPORTS1.1 reference databases must be built beforehand.  
