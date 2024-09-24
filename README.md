# miRNA.design.2024

This repository contains code for the design, analysis, and evaluation of a library of 3'UTRs containing miRNA target sites. [ADD LINK ONCE THE PAPER IS PUBLISHED]

## Table of Contents

1. [Library NGS Evaluation](#library-ngs-eval)
2. [Library Design](#library-design)
3. [Library Analysis](#library-analysis)

## Library NGS Eval

Evaluate the raw sequencing data to get counts and stabilities per construct.

### Requirements

- BWA-MEM
- PyDESeq2

### Usage

1. Download the data from GEO.
2. Set library names individually in each file (list at the top).
3. Extract the UMI from the reads and add it to the title:
   ```
   python 1_split_UTR_rem_UMIs.py
   ```
4. Generate reference files:
   - Use `2_ref_to_fasta.ipynb`
5. Index reference files:
   ```
   bwa index ./3UTR/references_3UTR.fasta
   ```
6. Align to the reference:
   ```
   nohup python3 -u 3_align_to_references.py > 3_align_to_references.txt &
   ```
7. Filter alignments:
   ```
   nohup python3 -u 4_filter_alignments.py > 4_filter_alignments.txt &
   ```
8. Sort and index:
   ```
   python 5_index_and_sort.py
   ```
9. UMI deduplication:
   ```
   python 6_deduplicate_umis.py
   ```
10. Count alignments:
    ```
    python 7_count_alignments.py
    ```
11. Process count data with PyDESEQ2:
    - Use `8_compute_lfc_DESeq2.ipynb`

**Note**: This process can also be run using a script in the folder. Set library names and generate index files before running.

## Library Design

Contains code for the design of library 2 and some analysis code for library 1 measurements.

### Requirements

- [NUPACK](https://docs.nupack.org/)

## Library Analysis

Code to generate all figures and analyses used in the paper.

### Requirements

- [NUPACK](https://docs.nupack.org/)

### Usage

To generate the paper figures:

1. Run the notebooks in order.
2. Some notebooks may need to be run multiple times or in a different order. This is indicated in the notebooks where necessary.