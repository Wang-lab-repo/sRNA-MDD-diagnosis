
# piRNA sRNA Pipeline


This pipeline describes how to quantify known piRNAs using Bowtie and RSEM, followed by extracting read-length distributions and merging expression values into a matrix.

---



## Step 1: Setup

```bash
cd /path/to/project
mkdir -p ./pirna/align ./pirna/pi_exp

clean_data=./clean
pirna_alin=./pirna/align

find $clean_data -name '*.fq' -exec basename {} \; | sed "s/\.fq//g" | sort > $pirna_alin/piRNA_dataset.txt
```

---

## Step 2: piRNA Mapping & Quantification

```bash
#!/bin/bash

base_dir=/path/to/project
clean_data=$base_dir/clean
pirna_alin=$base_dir/pirna/align
pi_exp=$base_dir/pirna/pi_exp
ref_dir=/path/to/reference/piRBase/piR_human

cat $pirna_alin/piRNA_dataset.txt | while read sample; do
  echo "Processing $sample..."

  bowtie -p 10 $ref_dir -v 1 -m 3 -q -S $clean_data/${sample}.fq     | samtools view -b -@ 4 > $pirna_alin/${sample}.bam

  rsem-calculate-expression --no-bam-output --alignments -p 10 --bam "$pirna_alin/${sample}.bam" "$ref_dir" "$pi_exp/${sample}.rsem"

  # Extract read-length and counts
  samtools view -F 4 "$pirna_alin/${sample}.bam" |     awk '{if($10 ~ /^[ATCGN]+$/) print $10}' |     awk 'BEGIN{OFS=","} {if(length($0) >= 15 && length($0) <= 45) {counts[$0]++; lengths[$0]=length($0)}} END {for (seq in counts) print seq, lengths[seq], counts[seq]}' |     sort -nr > "$pirna_alin/${sample}_read_count_length.csv"
done
```

---

## Step 3: Construct Expression Matrix

```r
library(dplyr)

folder_path <- "pirna/pi_exp"
files <- list.files(folder_path, pattern = "results$", full.names = TRUE)

output_folder <- file.path(folder_path, "../converted_counts")
dir.create(output_folder, showWarnings = FALSE)

for (file in files) {
  data <- read.table(file, header = TRUE, sep = "	")
  if (!("TPM" %in% colnames(data)) || !("expected_count" %in% colnames(data))) next
  total_expected <- sum(data$expected_count, na.rm = TRUE)
  data$calculated_count <- data$TPM * total_expected / 1e6
  output_name <- paste0(tools::file_path_sans_ext(basename(file)), "_calculated.tsv")
  output_path <- file.path(output_folder, output_name)
  write.table(data[, c("gene_id", "TPM", "expected_count", "calculated_count")],
              file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
}
message("✅ All files processed.")
```

---

## Step 4: Merge Count Matrix

```r
library(dplyr)
library(purrr)

input_folder <- "pirna/converted_counts"
files <- list.files(input_folder, pattern = "_calculated\.tsv$", full.names = TRUE)
count_list <- list()

for (file in files) {
  sample_name <- tools::file_path_sans_ext(basename(file)) %>% gsub("_calculated", "", .)
  data <- read.table(file, header = TRUE, sep = "\t")
  count_list[[sample_name]] <- data %>%
    select(gene_id, calculated_count) %>%
    rename(!!sample_name := calculated_count)
}

merged_counts <- reduce(count_list, full_join, by = "gene_id")
write.csv(merged_counts, file = file.path(input_folder, "merged_piRNA_counts.csv"),
          row.names = FALSE, quote = FALSE)
message("✅ Merged expression matrix saved.")
```

---

## Notes

- Ensure all reference paths are updated according to your cluster setup.  
- bowtie reference databases must be built beforehand.  

## Software Versions

| Tool            | Version     |
|-----------------|-------------|
| bowtie          | 1.0.0       |
| Trim Galore     | 1.2.28      |
| samtools        | 1.5         |