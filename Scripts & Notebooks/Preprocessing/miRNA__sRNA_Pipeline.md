# miRNA Quantification Pipeline

This document outlines the standalone pipeline used to quantify known miRNAs from clean small RNA sequencing data using [miRDeep2](https://github.com/rajewsky-lab/mirdeep2).  

---

## Step 1: Setup

```bash
cd /path/to/project
mkdir -p mirna

# Full path list
find ./clean -name '*_trimmed.fq' > mirna_name.list

# Short sample names
find ./clean -name '*_trimmed.fq' -exec basename {} \; | sed 's/_trimmed\.fq$//' > ./mirna/short_mirna_name.list

```

---

## Step 2: miRNA Mapping & Quantification

```bash
#!/bin/bash

export PATH=/path/to/anaconda3/bin:$PATH

# Define paths
base_dir=/path/to/project
clean_data=$base_dir/clean
mirna_data=$base_dir/mirna
genome_dir=/path/to/reference/Homo_sapiens
mirdeep_dir=/path/to/mirdeep2/bin

# Run in loop for all samples
cat mirna/short_mirna_name.list | while read sample; do
    echo "Processing $sample..."

    # Mapping reads to genome using mapper.pl
    perl $mirdeep_dir/mapper.pl $clean_data/${sample}.fq \
        -e -h -j -k AGATCGGAAGAGC -l 18 -m \
        -p $genome_dir/genome/genome \
        -s $mirna_data/${sample}_mapped.fa \
        -t $mirna_data/${sample}_vs_genome.arf -v

    # Quantifying known miRNAs using quantifier.pl
    perl $mirdeep_dir/quantifier.pl \
        -p $genome_dir/miRBase/hairpin.fa \
        -m $genome_dir/miRBase/mature.fa \
        -r $mirna_data/${sample}_mapped.fa \
        -t hsa -y ${sample}_mirna
done
```

---

## Step 3: Construct Expression Matrix


```bash
cd /path/to/project
mkdir -p mirna_merge/input

find ./mirna -type f -name "*.csv" > mirna_name.list
find ./mirna -type f -name "*.csv" -exec basename {} \; > mirna_short_name.list

cat mirna_name.list | while read id; do
  if [[ -f $id ]]; then
    file_name=$(basename "$id" | cut -d '_' -f 5)
    echo "<=== $file_name ===>"

    awk 'BEGIN {OFS="\t"} {printf "%s_%s\t%s\n", $1, $3, $2}' "$id" | sed '1d' > "./mirna_merge/input/${file_name}.csv"

    printf "mirna_name_precursor\t%s\n" "$file_name" > "./mirna_merge/input/${file_name}_header.csv"
    cat "./mirna_merge/input/${file_name}.csv" >> "./mirna_merge/input/${file_name}_header.csv"
  else
    echo "Warning: File $id not found."
  fi
done
```

---

## Step 4: Merge Into Final Matrix

```bash
cd mirna_merge
python mirna_merge.py
```
---

## Reference Databases

Genome: Human hg19  

miRNA Database:  

Hairpin and mature sequences from miRBase  

Downloaded FASTA: hairpin.fa and mature.fa  

---

## Notes

- Ensure all reference paths are updated according to your cluster setup.  
- miRDeep2 reference databases must be built beforehand.  