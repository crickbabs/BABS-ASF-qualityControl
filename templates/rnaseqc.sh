#!/bin/sh

# the conda wrapper does not allow us to pass all the arguments provided by
# RNA-SeQC, we just get the path of the jar file and run our own command

# can fail otherwise
export LC_ALL=en_US.utf8
export LANG=en_US.utf8

# i suspect the conda wrapper follows links...
local_bam=\$(basename $bam .bam).copy.bam
local_bai=\$(basename $bai .bam.bai).copy.bam.bai
cp -Lv $bam ./\$local_bam
cp -Lv $bai ./\$local_bai

rnaseqc \
	--reference $fasta \
	--annotation-gtf $gtf \
	--sample-name  $basename \
	--output-dir $metrics_filename \
	--dry-run \
	\$local_bam \
	> command.txt

jar_path=\$(sed "s/.*'\\(\\/.*\\.jar\\)'.*/\\\\1/" command.txt)

java -Xmx10g -jar \$jar_path \
	-d 1000000 \
	-rRNA $rrna \
	-r $fasta \
	-t $gtf \
	-o $metrics_filename \
	-gatkFlags '-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY' \
	$single_end_arg \
	-s "$basename|\$local_bam|$basename"

# clean
rm -v \$local_bam
rm -v \$local_bai

