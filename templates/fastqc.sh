#!/bin/sh

mkdir $outdir

# to avoid multiple the multiqc general statistics table
cp -Lv $fastq $new_fastq

fastqc --outdir $outdir $new_fastq
rm -v $new_fastq

