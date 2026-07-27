# Small RNA Preprocessing Pipeline

This project provides a two-step pipeline for quality control and trimming of small RNA sequencing data.

Licensed under the **Apache License, Version 2.0**.

---

## Step 1: QC and Trimming

```bash
#!/bin/bash

base_dir="/path/to/project"
cd "${base_dir}"
mkdir -p sequence trim fastqc clean

raw_data="${base_dir}/sequence"
find "${raw_data}" -type f -name "*.gz" > "${base_dir}/name.list"
find "${raw_data}" -type f -name "*.gz" -exec basename {} \; > "${base_dir}/short_name.list"

while read -r id; do
    echo "<=== ${id} ===>"
    md5sum "${id}"
done < "${base_dir}/name.list" > md5.log

export PATH="/path/to/parallel/bin:/path/to/anaconda3/bin:/path/to/TrimGalore:$PATH"

while read -r id; do
    filename=$(basename "${id}")
    prefix="${filename%%.*}"

    fastqc -t 8 -o "${base_dir}/fastqc" "${id}"

    trim_galore \
        --quality 25 \
        --phred33 \
        --fastqc \
        --length 18 \
        --stringency 3 \
        -o "${base_dir}/trim" "${id}" > "${base_dir}/trim/${prefix}_trim.log" 2>&1
done < "${base_dir}/name.list"
```

---

## Step 2: Advanced Filtering by Base Quality

```bash
python filter_by_quality.py input.fq.gz ./clean/filtered_output.fq
```
---

## Software Versions

| Tool            | Version     |
|-----------------|-------------|
| FastQC          | 0.12.1     |
| Trim Galore     | 0.6.7      |
| Python          | 3.11       |