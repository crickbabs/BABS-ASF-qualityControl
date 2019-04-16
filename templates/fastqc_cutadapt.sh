#!/bin/sh

mkdir $outdir
fastqc --outdir $outdir $fastq

