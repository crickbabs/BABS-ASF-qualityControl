#!/bin/sh

# this the command inside a conda environment

picard \
	-Xmx10g \
	-Djava.io.tmpdir=$tmp_dirname \
	AddOrReplaceReadGroups \
		VALIDATION_STRINGENCY=SILENT \
		INPUT=$bam \
		OUTPUT=$filename \
		RGID=$name \
		RGLB=$name \
		RGPU=$name \
		RGSM=$name \
		RGCN=TheFrancisCrickInsitute \
		RGPL=Illumina

