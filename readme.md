# GO Enrichment – Promoter Region Extraction

## Project Overview

This project focuses on extracting promoter regions from human gene annotation data using strand-aware genomic coordinates.
The promoter regions generated here represent **500 bp upstream of the Transcription Start Site (TSS)** and can be used for:

* GO enrichment analysis
* Motif discovery
* Regulatory sequence analysis
* Functional genomics studies

The workflow was performed using `bedtools`, `samtools`, and simple command-line processing.

---

# Files Used

| File Name              | Purpose                                         |
| ---------------------- | ----------------------------------------------- |
| `Genes_TSS_final.bed`  | Final cleaned BED file containing TSS positions |
| `promoters_500_bp.bed` | Generated 500 bp promoter regions               |

---

# Software and Environment

Environment created using Conda/Mamba:

```bash
mamba create -n go_enrichment python=3.12
mamba activate go_enrichment
```

Packages installed:

```bash
mamba install -c bioconda bedtools samtools emboss
```

---

# Input Data

## Reference Genome

Human genome assembly (hg38) downloaded from UCSC:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
```

Genome FASTA file:

```text
hg38.fa
```

---

## Annotation File

Gene annotation file used:

```text
human_gene_annotation.tsv.gz
```

---

# Step 1 — Prepare Genome Index

Index the reference genome:

```bash
samtools faidx hg38.fa
```

Create chromosome size file:

```bash
cut -f1,2 hg38.fa.fai > hg38.genome
```

---

# Step 2 — Create TSS BED File

The annotation file was converted into BED format using `awk`.

Command used:

```bash
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"}
NR>1{

    if($8==-1 || $8=="" || $7=="")
        next

    chrom="chr"$5

    if(chrom=="chrMT")
        chrom="chrM"

    tss=$8 + 0

    strand="+"
    if($6==-1)
        strand="-"

    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand
}' > genes_tss_clean.bed
```

Remove chromosomes not present in the genome file:

```bash
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > Genes_TSS_final.bed
```

---

# BED File Format

The generated BED file contains:

| Column | Information                     |
| ------ | ------------------------------- |
| 1      | Chromosome                      |
| 2      | Start coordinate                |
| 3      | End coordinate                  |
| 4      | Gene identifier                 |
| 5      | Placeholder (`.`)               |
| 6      | Strand information (`+` or `-`) |

---

# Step 3 — Generate Promoter Regions

Promoter regions were generated using `bedtools slop`.

Command:

```bash
bedtools slop \
-i Genes_TSS_final.bed \
-g hg38.genome \
-l 500 \
-r 0 \
-s \
> promoters_500_bp.bed
```

---

# Strand-Aware Promoter Extraction

The `-s` option in `bedtools slop` makes the promoter extraction strand-aware.

This means:

* For `+` strand genes → promoter extends toward smaller coordinates
* For `-` strand genes → promoter extends toward larger coordinates

Example:

```text
chrM    4400    4901    chrM@4400-4401|MT-TQ    .    -
```

Here, the promoter region was correctly extended according to the negative strand direction.

---

# Final Output

Final promoter file:

```text
promoters_500_bp.bed
```

This file contains strand-aware promoter regions that can be used for downstream GO enrichment and regulatory analysis workflows.

