#!/bin/sh

seqtk sample -s1903 ${fastq} ${n} | gzip -c > ${filename}

