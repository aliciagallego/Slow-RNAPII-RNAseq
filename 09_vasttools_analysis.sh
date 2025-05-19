#!/bin/bash

# This script computes vast-tools analysis following author's instructions in: 
# https://github.com/vastgroup/vast-tools

# set working directory
cd /path/Splicing/vasttools/

# number of threads
THREADS=$(( $(nproc) - 2 ))

# 1. vast-tools align
# -------------------
# This function computes PSI values for each provided sample

# fastq sample names (rename accordingly) 
for SAMPLE in sample1 sample2 sample3 sample4; do
	INPUT="/input_path/"$SAMPLE
	vast-tools align $INPUT"/"$SAMPLE"_1.fq.gz" $INPUT"/"$SAMPLE"_2.fq.gz" -sp mm10 --expr -c $THREADS
done 

# 2. vast-tools merge
# -------------------
# This function can be used to increase read coverage. It requires a configuration file (config_file) to specify which replicates should be merged
vast-tools merge --groups config_file --move_to_PARTS --sp mm10

# 3. vast-tools combine
# ---------------------
# This function generates one final table named as INCLUSION_LEVELS_FULL-mm10-(number_of_samples).tab
vast-tools combine -sp mm10 --cores $THREADS

# 4. vast-tools compare 
# ---------------------
# This function computes dPSIs between pairs of samples

# --min_dPSI: minimum absolute value of ΔPSI (default >15).
# --min_range: Minimum distance between the ranges of both groups (default >5). 
# --outFile file: output file name (default based on option parameters).
# -a/--samplesA sA1,sA2    Required, 1:n sample names or column_# separated by , (mandatory) (usually WT)
# -b/--samplesB sB1,sB2    Required, 1:n sample names or column_# separated by , (mandatory) (usually KO)
# --only_samples: to plot only the compared samples (by default, all samples in the INCLUSION file are plotted).

# define INCLUSION_TABLE.tab obtained in vast-tools combine
INCLUSION_TABLE="/media/cc/C/Joana/totalRNA/totalRNA_output/12_Splicing/vasttools/vast_out/INCLUSION_LEVELS_FULL-mm10-*.tab"

# define samples -a and -b accordingly 
# example: naive vs primed comparisons
vast-tools compare $INCLUSION_TABLE -a WTd0  -b WTd1 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10
vast-tools compare $INCLUSION_TABLE -a WTd0  -b WTd2 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10
vast-tools compare $INCLUSION_TABLE -a WTd0  -b WTd7 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10
vast-tools compare $INCLUSION_TABLE -a mutd0  -b mutd1 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10
vast-tools compare $INCLUSION_TABLE -a mutd0  -b mutd2 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10
vast-tools compare $INCLUSION_TABLE -a mutd0  -b mutd7 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10

# example: mut vs WT comparison
vast-tools compare $INCLUSION_TABLE -a WTd0  -b mutd0 --min_dPSI 15 --min_range 5 --print_dPSI --only_samples --sp mm10
