#!/bin/sh

module load Nextflow/0.32.0

# work directory
WORK=/camp/stp/babs/scratch/bahn
WORK=$WORK/qcpipeline
WORK=$WORK/work

nextflow run main.nf \
	-lib qcpipeline.jar \
	-params-file params.yml \
	-work-dir $WORK \
	-with-timeline nextflow/timeline.html \
	-with-report nextflow/report.html \
	-with-dag nextflow/dag.dot \
	-resume

