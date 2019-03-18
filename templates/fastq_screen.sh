#!/bin/sh

fastq_screen \
	--force \
	--outdir ./ \
	--subset 200000 \
	--conf ${fastq_screen_conf}  \
	--threads ${task.cpus} \
	--aligner bowtie2 \
	${fastq}

