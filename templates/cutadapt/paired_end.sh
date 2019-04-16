#!/bin/sh

# to avoid multiple the multiqc general statistics table
cp -v $fastq1 $new_fastq1
cp -v $fastq2 $new_fastq2

cutadapt \
	-a AGATCGGAAGAGC \
	-A AGATCGGAAGAGC \
	-o $output1 \
	-p $output2 \
	-e 0.1 \
	-q 10 \
	-m 25 \
	-O 1 \
	$new_fastq1 $new_fastq2 > $logfile

rm -v $new_fastq1
rm -v $new_fastq2

