#!/bin/bash

# fastq sample names (rename accordingly) 
for SAMPLE in sample1 sample2 sample3 sample4; do

# paired-end files
for NUM in 1 2; do

	INPUT="/input_path/"$SAMPLE"/"$SAMPLE"_"$NUM".fq.gz"
	OUTPUT="/output_path/"
	fastqc -o $OUTPUT $INPUT
done
done 
