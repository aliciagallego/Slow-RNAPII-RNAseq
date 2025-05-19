#!/bin/bash

# bam file names (rename accordingly) 
for SAMPLE in sample1 sample2 sample3 sample4; do
	
	# path of bam files
        BAM="/input_path/"$SAMPLE"Aligned.sortedByCoord.out.bam" 
	
	# output path
       	OUTFW="/output_path/"$SAMPLE"_FW.bw"
       	OUTRV="/output_path/"$SAMPLE"_RV.bw"
	
	# do bam index
	samtools index $BAM
       	
       	# do bigwig files
       	bamCoverage -b $BAM -bs 1 --normalizeUsing CPM --filterRNAstrand forward -o $OUTFW
       	bamCoverage -b $BAM -bs 1 --normalizeUsing CPM --filterRNAstrand reverse --scaleFactor -1 -o $OUTRV      
done
