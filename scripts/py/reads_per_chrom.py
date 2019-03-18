#!/usr/bin/python

import sys, re
import json
import pandas as pd

idxstats_filepath = sys.argv[1]

# load
hdr = ["sequence_name", "sequence_length", "mapped_reads", "unmapped_reads"]
df = pd.read_csv(idxstats_filepath, sep="\t", names=hdr)

# filter
df = df.loc[
		df["sequence_name"].\
				apply( lambda x: bool(re.match("(^\d+$)|^X$|^Y$|^MT$", x)) )
		]

# compute proportions
df["proportion"] = df["mapped_reads"].\
		apply( lambda x: x / sum(df["mapped_reads"]) )
df["expected_proportion"] = df["sequence_length"].\
	apply( lambda x: x / sum(df["sequence_length"]) )
df["gap"] = abs( df["proportion"] - df["expected_proportion"] )
df["diff"] = df["proportion"] - df["expected_proportion"]

# compute counts
df["expected_count"] = df["expected_proportion"].\
		apply( lambda x: x * sum(df["mapped_reads"]) )
df["gap_count"] = abs( df["mapped_reads"] - df["expected_count"] )
df["diff_count"] = df["mapped_reads"] - df["expected_count"]

# sort (very dirty :/)
df["sequence_name"] = df["sequence_name"].\
		apply(lambda x: re.sub("^(\d)$", "0\\1", x))
df = df.sort_values("sequence_name")
df["sequence_name"] = df["sequence_name"].\
		apply(lambda x: re.sub("^0", "", x))

# export
print(json.dumps( df.to_dict("list"), indent=3 ))

