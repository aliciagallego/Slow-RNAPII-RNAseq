#!/bin/bash

# STAR alignment - parameters for repeats

# genome index path
INDEX="/path/STAR_index_mm10/"

# number of threads
THREADS=$(( $(nproc) - 2 ))

# fastq sample names (rename accordingly) 
for SAMPLE in sample1 sample2 sample3 sample4; do

	FASTQ="/input_path/*${SAMPLE}*"
	STAR --runMode alignReads \
	     --genomeDir $INDEX \
	     --readFilesIn $FASTQ/*_1.fq.gz $FASTQ/*_2.fq.gz \
	     --readFilesCommand zcat \
	     --outSAMtype BAM SortedByCoordinate \
	     --outFileNamePrefix $SAMPLE \
	     --outFilterType Normal \
	     --alignEndsType Local \
	     --winAnchorMultimapNmax 5000 \
	     --seedPerWindowNmax 1000 \
	     --alignTranscriptsPerReadNmax 100000 \
	     --seedNoneLociPerWindow 100 \
	     --alignWindowsPerReadNmax 20000 \
	     --alignTranscriptsPerWindowNmax 1000 \
	     --outFilterMultimapNmax 100000 \
	     --outSAMmultNmax 1
done
