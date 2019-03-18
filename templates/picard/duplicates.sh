#!/bin/sh

# this the command inside a cond environment

picard \
	-Xmx10g \
	-Djava.io.tmpdir=$tmp_dirname \
	MarkDuplicates \
		VALIDATION_STRINGENCY=SILENT \
		INPUT=$bam \
		OUTPUT=$filename \
		METRICS_FILE=$metrics_filename \
		ASSUME_SORTED=true \
		REMOVE_DUPLICATES=false \
		TMP_DIR=$tmp_dirname

