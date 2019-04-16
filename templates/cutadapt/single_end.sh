#!/bin/sh

# to avoid multiple the multiqc general statistics table
cp -Lv $fastq $new_fastq

cutadapt \
	-a AGATCGGAAGAGC \
	-o $output \
	-e 0.1 \
	-q 10 \
	-m 25 \
	-O 1 \
	$new_fastq > $logfile

rm -v $new_fastq

