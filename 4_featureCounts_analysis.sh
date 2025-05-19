#!/bin/bash

# Description of parameters
# -------------------------
# -a: genome annotation file (example_genome_annotation.gtf): use the same GTF file that was used in the STAR alignment
# -T: number of the threads. 1 by default.
# -p: specify that input data contain paired-end reads. 
# -s: strand-specific read counting (library used in the sequencing protocol). 0 = unstranded reads, 1 = stranded reads, 2 = reversely stranded reads.
# -M: multi-mapping reads/fragments will be counted.
# -o: name of the output file.
# -B: only fragments that have both ends successfully aligned will be considered for summarization. This option is only applicable for paired-end reads.
	
# input path to bam files
BAMS="/input_path/*.sortedByCoord.out.bam"
	
# output path
OUT="/output_path/featureCounts_file.txt"
	
# genome index path
GTF="/path/STAR_index_mm10/mm10_genome_annotation.gtf"
	
# do featureCounts
featureCounts -t transcript -a $GTF $BAMS -T 8 -p -s 2 -B -o $OUT
