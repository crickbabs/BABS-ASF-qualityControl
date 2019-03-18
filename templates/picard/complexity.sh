#!/bin/sh

picard \
	-Xmx10g \
	-Djava.io.tmpdir=$tmp_dirname \
	EstimateLibraryComplexity \
		VALIDATION_STRINGENCY=SILENT \
		INPUT=$bam \
		OUTPUT=$metrics_filename \
		TMP_DIR=$tmp_dirname

