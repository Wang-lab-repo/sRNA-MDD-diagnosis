# sRNA Target Prediction Pipeline

This pipeline describes how to predict sRNA targets from clean sRNA-seq fasta file using standard tools like miRanda, RNAhybrid, and post-processing in R/Python.

## 1. Input Preparation

- Cleaned FASTA files of small RNA.  

```bash
cat mirna_seq.fa rsrna.fa tsrna.fa pirna_seq.fa > seq.fa

# >piR-hsa-114184
# GGCTGGTCAGATGGTAGTGGGTTATCAGAACTT
# >piR-hsa-1249873
# CGGACCAGGGGAATCCGACTGTCTAATTAAAAC
# >piR-hsa-138201
# GGCTGGTCCGATGGTAGTGGGTTCTCAGAACTT
```
- Reference in FASTA format.  

## 2. Target Prediction with miRanda

### Run miRanda:

```bash
miranda seq.fa UTR_clean.fasta -out miranda_output.txt -quiet
grep '>>' miranda_output.txt > miranda_result.txt
```

### Parse miRanda Output:

Use `run_miranda_parser.py` to extract results into a clean table.

```bash
python run_miranda_parser.py
```

---

## 3. Target Prediction with RNAhybrid

### Run RNAhybrid:

```bash
RNAhybrid -g jpg -b 1 -e -20 -f 8,12 -u 1 -v 1 -s 3utr_human -t UTR_clean.fasta -q seq.fa > rnahybrid_output.txt
```

### Parse RNAhybrid Output:

```bash
python parse_rnahybrid.py
```

---

## 4. Intersect Targets from miRanda & RNAhybrid

To increase prediction confidence, take intersection of miRanda and RNAhybrid results.

```bash
python intersect_targets.py
```

---

## Software Versions

| Tool            | Version     |
|-----------------|-------------|
| miranda         | v3.3a       |
| RNAhybrid    | v2.1.2      |


