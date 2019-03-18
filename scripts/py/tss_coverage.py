#!/usr/bin/python

# coding: utf-8

# nourdine.bah@crick.ac.uk

# this code has been taken there:
# http://htseq.readthedocs.io/en/master/tss.html
# thanks to them

from os.path import join
from os.path import exists
from os.path import basename
import argparse
import HTSeq
import itertools
import numpy as np
import json


def compute_tss_coverage(
		bam_filepath,
		gtf_filepath,
		half_winwidth=3000,
		fragmentsize=200):

	# htseq object creation
	assert exists(gtf_filepath), "The provided GTF file path does not exist."
	assert exists(bam_filepath), "The provided BAM file path does not exist."
	gtf = HTSeq.GFF_Reader(gtf_filepath)
	bam = HTSeq.BAM_Reader(bam_filepath)
	
	# all the TSSs
	tsspos = set()
	#for feat in itertools.islice(gtf, 1000): # for debug
	for feat in gtf:
		if feat.type == "exon" and feat.attr["exon_number"] == '1':
			tsspos.add(feat.iv.start_d_as_pos)
	
	# nicer like that
	tsspos = sorted(list(tsspos), key=lambda x: (x.chrom, x.start))

	support = np.arange( -half_winwidth, half_winwidth )
	profile = np.zeros( 2 * half_winwidth, dtype='i')

	for p in tsspos:

		window = HTSeq.GenomicInterval(
				p.chrom,
				max(p.pos - half_winwidth - fragmentsize, 0),
				p.pos + half_winwidth + fragmentsize,
				".")

		for a in bam[window]:
			a.iv.length = fragmentsize

			if p.strand == "+":
				start_in_window = a.iv.start - p.pos + half_winwidth
				end_in_window = a.iv.end - p.pos + half_winwidth
			else:
				start_in_window = p.pos + half_winwidth - a.iv.end
				end_in_window = p.pos + half_winwidth - a.iv.start

			start_in_window = max(start_in_window, 0)
			end_in_window = min(end_in_window, 2*half_winwidth)

			if start_in_window >= 2*half_winwidth or end_in_window < 0:
				continue
			else:
				profile[start_in_window:end_in_window] += 1
	
	return ( support.tolist(), profile.tolist() )
	

if __name__=="__main__":
	
	description = "Returns a TSS coverage as a json."
	ap = argparse.ArgumentParser(description=description)

	# a bam file, an annotation, a window size and a fragment size
	ap.add_argument("-b", "--bam", help="the BAM file", required=True)
	ap.add_argument("-g", "--gtf", help="the GTF file", required=True)
	ap.add_argument("-w", "--win",
			help="the half window width in bp",
			required=False,
			default="3000")
	ap.add_argument("-s", "--size",
			help="the fragment size in bp",
			required=False,
			default="200")
	ap.add_argument("-n", "--name", help="the label", required=False)

	args = ap.parse_args()

	sup, cov = compute_tss_coverage(
			args.bam,
			args.gtf,
			int(args.win),
			int(args.size)
			)

	# need to add a tsv export so the user can have the choice
	name = args.name if args.name else basename(args.bam)
	export = { "name": name, "support": sup, "coverage": cov }
	#print( json.dumps(export, indent=3) )
	print( json.dumps(export) )


