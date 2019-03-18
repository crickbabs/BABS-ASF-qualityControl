#!/bin/sh

STAR \
	--runMode alignReads \
	--runThreadN $cpus \
	--outSAMtype BAM Unsorted SortedByCoordinate \
	--outSAMattrRGline \
		ID:$basename \
		LB:$basename \
		PU:$basename \
		SM:$basename \
		CN:TheFrancisCrickInsitute \
		PL:Illumina \
		DS:BABS_QC_pipeline \
	--outFileNamePrefix $basename. \
	--genomeDir $genome \
	--readFilesCommand zcat \
	--limitBAMsortRAM 6000000000 \
	--readFilesIn $fastq

