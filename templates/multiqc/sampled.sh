#!/bin/sh

export LANG=en_US.utf8
export LC_ALL=en_US.utf8

multiqc \
	--force \
	--data-format $multiqc_data_format \
	--template $multiqc_template \
	--config $multiqc_conf \
	.

