# miRNA.design.2025

This repository contains code for the design, analysis, and evaluation of a library of 3'UTRs and 5'UTRs containing miRNA target sites.
It also contains easily adaptable example code for how to unify microRNA data from two sources, filter it for crosstalk and other confounders,
then use it for prediction and design of microRNA target sites with desired expression properties.

## Contents
### example_processing_and_design
Example code for how to process microRNA data to merge data sources, remove crosstalk and other confounders, and make predictions and designs.
The three notebooks should be run in order. Input data is provided. This is the only part necessary for most users.

### libraries_ngs_data_processing
NGS analysis code to map NGS counts to UMIs per construct

### library_design
Code for used to design Library 2. Code for other libraries is analogous.

### library_ratiometric_eval
Analysis code for the ratiometric stability data

### library_actd_and_polysome_analysis
Analysis code for Actinomycin D and Polysome Profiling
Analysis code for the ratiometric stability data

## Requirements and Installation
All code except the NGS processing code only requires standard python packages. As such, it should run on essentially any machine capable of executing python code.
No special installation is necessary.

NGS data analysis requires:
BWA: https://github.com/lh3/bwa
UMI-tools: https://umi-tools.readthedocs.io/en/latest/INSTALL.html

## Usage
1. Generally, run the notebooks in a subfolder in the indicated order (1_xx.ipynb, 2_xx.ipynb ...).
2. For the results of the publication itself, some notebooks may need to be run multiple times or in a different order. This is indicated in the notebooks where necessary.