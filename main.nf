#!/usr/bin/env nextflow

/*
 * christopher.barrington chez crick.ac.uk
 * gavin.kelly chez crick.ac.uk
 * harshil.patel chez crick.ac.uk
 * nourdine.bah chez crick.ac.uk
 * philip.east chez crick.ac.uk
 * richard.mitter chez crick.ac.uk
 */

// groovy modules
import groovy.io.FileType
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

// the java modules
import java.nio.file.Paths

// custom module
import qcpipeline.Project
import qcpipeline.Sample
import qcpipeline.Fastq

publish_mode = "symlink"
publish_overwrite = true

// all the fastq files of the project
fastq_files = []
params.directories.each{ path ->
	def dir = new File(path)
	dir.eachFileMatch(FileType.FILES, ~/.*\.fastq\.gz/) { f ->
		Fastq fastq = new Fastq(f)
		fastq.setSpecies(params.species)
		fastq.setSequencing_type(params.type.toLowerCase().replaceAll(" |-", ""))
		fastq_files.add(fastq)
	}
}

/*|||||||||||||||||||||*/
Channel.from(fastq_files)
/*|||||||||||||||||||||*/
	//// test
	//.map{
	//	[ it.getFile().getBaseName().replaceAll("_R\\d_001.fastq", ""), it ]
	//}
	//.groupTuple()
	//.take(3)
	//.map{ it[1] }
	//.flatten()
	//// test
	.into{ fastqs; to_read_length }

/*||||||||||*/
to_read_length
/*||||||||||*/
	.map{ [ it.getBasename() , it.getFile() ] }
	.set{ to_read_length }

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ---------------------------- PARAMETERS --------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

// choose the species
genomes = [:]
genomes["homo sapiens"] = ["version":"GRCh38", "release":"86"]
genomes["mus musculus"] = ["version":"GRCm38", "release":"86"]

// configuration
fastq_screen_conf = absolute_path( "conf/fastq_screen.conf" )
multiqc_conf_file = absolute_path( "conf/multiqc.yml" )

// script locations
py_dir = absolute_path("scripts/py")
cpp_dir = absolute_path("scripts/cpp")
qc_genes_dir = absolute_path("scripts/qc_genes_list/txt")

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------- VARIOUS FUNCTIONS ----------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*--------------------------- */
def get_file(species, name) {
/*--------------------------- */

	version = genomes[ species.toLowerCase() ]["version"]
	release = genomes[ species.toLowerCase() ]["release"]

	genome = params.genomes[version][release]

	return genome[name]
}

/* ----------------------------------------------- */
def get_star_index(species, rough_read_length) {
/* ----------------------------------------------- */

	version = genomes[ species.toLowerCase() ]["version"]
	release = genomes[ species.toLowerCase() ]["release"]

	genome = params.genomes[version][release]

	return genome.star[ rough_read_length + "bp" ].index
}

/* ----------------------------------------------- */
def get_rsem_star_index(species, rough_read_length) {
/* ----------------------------------------------- */

	version = genomes[ species.toLowerCase() ]["version"]
	release = genomes[ species.toLowerCase() ]["release"]

	genome = params.genomes[version][release]

	return genome.rsem.rsem_star[ rough_read_length + "bp" ].index
}

/* --------------------------- */
def flatten_paired_end(channel) {
/* --------------------------- */

	// the .transpose() nextflow operator can be used here instead

	def single_end = channel[2]

	if (single_end) {

		return channel

	} else {

		return [ *channel[0..2], channel[3][0], *channel[0..2], channel[3][1] ]

	}
}

/* ------------------------ */
def sort_files(path1, path2) {
/* ------------------------ */

	def f1 = Paths.get( path1.toString() ).getFileName()
	def f2 = Paths.get( path2.toString() ).getFileName()

	if ( f1 == f2 ) {
		return 0
	} else if ( f1 > f2 ) {
		return 1
	} else {
		return -1
	}
}

/* ------------------- */
def absolute_path(path) {
/* ------------------- */

	def f = new File(path)

	if ( ! f.isAbsolute() ) {
		return Paths.get(workflow.projectDir.toString(), path).toString()
	} else {
		return path
	}
}

/* -------------- */
def get_lims(name) {
/* -------------- */

	return name.replaceAll("_.*", "")
}

/* -------------------------- */
def get_experiment_info(fastq) {
/* -------------------------- */

	// fastq info
	def species = fastq.getSpecies()
	def read_length = fastq.getRead_length()
	def rough_read_length = fastq.getRough_read_length()
	def type = fastq.getSequencing_type()
	def machine = fastq.getMachine().getName()
	def run = fastq.getRun().getNumber()
	def lane = fastq.getLane().getNumber()


	// annotation
	def star_index = get_star_index( species , rough_read_length )
	def rsem_star_index = get_star_index( species , rough_read_length )
	def fasta = get_file( species , "fasta" )
	def bed = get_file( species , "bed" )
	def refflat = get_file( species , "refflat" )
	def rrna_list = get_file( species , "rrna_list" )
	def rrna_interval_list = get_file( species , "rrna_interval_list" )
	def gtf = get_file( species , "gtf" )
	def rnaseqc_gtf = get_file( species , "rnaseqc_gtf" )

	def info = [
		[ "Species" , species ],
		[ "Read Length" , read_length ],
		[ "Rough Read Length" , rough_read_length ],
		[ "Sequencing Type" , type ],
		[ "Machine" , machine ],
		[ "Run" , run ],
		[ "Lane" , lane ],
		[ "STAR Index" , star_index ],
		[ "RSEM-STAR Index" , rsem_star_index ],
		[ "FASTA" , fasta ],
		[ "GTF" , gtf ],
		[ "BED" , bed ],
		[ "REFFLAT" , refflat ],
		[ "RNASeqC GTF" , rnaseqc_gtf ],
		[ "Ribosomal RNAs" , rrna_list ],
		[ "Ribosomal RNAs Intervals" , rrna_interval_list ]
	]

	return info
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ---------------------------- READ LENGTH -------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

///////////////////
process read_length {
///////////////////

	tag { basename }

	input:
		set val(basename), val(fastq) from to_read_length

	output:
		set val(basename), stdout into lengths

	script:
		template "read_length.sh"
}

/*||*/
fastqs
/*||*/
	.map{ [ it.getBasename() , it ] }
	.combine(lengths)
	.filter{ it[0] == it[2] }
	.map{ it[1].addReadLength( it[3].toInteger() ) }
	.set{ fastqs }

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// --------------------------- CREATE PROJECT ------------------------------ //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

// all the fastq files with their read length
fastq_list =
	fastqs
		.toSortedList()
		.get()
		.toList()

// finally create the project
project = new Project(fastq_list)
project_outdir = project.getId()

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------- FASTQC ---------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

//////////////
process fastqc {
//////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(fastq) from Channel.from( project.getFastqC() )

	output:
		set \
			val(flowcell),
			val(name),
			file("*.fastqc") into fastqc
	
	script:
		basename = flowcell + "_" + name
		outdir = basename + ".fastqc"
		template "fastqc.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------ CUTADAPT --------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////
process cutadapt {
////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(fastq) from Channel.from( project.getCutadapt() )
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file('*.log') into cutadapt_log
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file('*.cutadapt.fastq.gz') \
			into \
				cutadapt,
				cutadapt_fastq_screen,
				cutadapt_align,
				cutadapt_multiqc

	script:

		basename = flowcell + "_" + name

		if (single_end) {

			fastq = fastq[0]
			output = basename + ".cutadapt.fastq.gz"
			logfile = basename + ".log"

			template "cutadapt/single_end.sh"

		} else {

			fastq1 = fastq[0]
			fastq2 = fastq[1]
			output1 = basename + "_1.cutadapt.fastq.gz"
			output2 = basename + "_2.cutadapt.fastq.gz"
			logfile = basename + ".log"

			template "cutadapt/paired_end.sh"
		}
}

/*|| separate the two fastq files if it's paired-end and add a tag ||*/
/*|| corresponding to the read number                              ||*/
/*||||*/
cutadapt
/*||||*/
	.map{ flatten_paired_end(it) } // .transpose() can be used here instead
	.flatten()
	.collate(4)
	.map{[
		*it[0..2],
		it[3] ==~ /.*_1.cutadapt.fastq.gz$/ ? "1" : "2",
		it[3]
	]}
	.into{ cutadapt_fastqc; cutadapt_strandedness }

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// -------------------------- FASTQC CUTADAPT ------------------------------ //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

///////////////////////
process fastqc_cutadapt {
///////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(read_num),
			file(fastq) from cutadapt_fastqc

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(read_num),
			file("*fastqc*") into fastqc_cutadapt
	
	script:
	
		basename = flowcell + "_" + name + "_R" + read_num
		basename = "cutadapt_" + basename
		outdir = basename + ".fastqc"

		template "fastqc.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// -------------------------- FASTQ SCREEN --------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////////
process fastq_screen {
////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(fastq) from cutadapt_fastq_screen.transpose()

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.html") into fscreen_html
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.txt") into fscreen_txt
	
	script:

		basename = flowcell + "_" + name

		template "fastq_screen.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// -------------------- STRANDEDNESS SEQUENCING TYPE  ---------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||||||||||||||||||*/
rnaseq_cutadapt_strandedness =
	Channel.create()

/*||||||||||||||||||||||||||||*/
non_rnaseq_cutadapt_strandedness =
	Channel.create()

/*|||||||||||||||||*/
cutadapt_strandedness
	.choice(rnaseq_cutadapt_strandedness, non_rnaseq_cutadapt_strandedness) {
		it = project.getSequencingType( get_lims(it[1]) ) == "rnaseq" ? 0 : 1
	}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ----------------------------- SAMPLING ---------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// ------------------------------ SEQTK ------------------------------------ //
///////////////////////////////////////////////////////////////////////////////

////////////////
process sampling {
////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(read_num),
			file(fastq) from rnaseq_cutadapt_strandedness

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(read_num),
			file("*.sampled.fastq.gz") into sampling
	
	script:

		n = 1000
		basename = flowcell + "_" + name + "_R" + read_num
		filename =  basename + ".sampled.fastq.gz"

		template "seqtk/sample.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------- STAR ------------------------------------ //
///////////////////////////////////////////////////////////////////////////////

/*||||*/
sampling
/*||||*/
	.groupTuple(by:[0,1,2])
	.map{ [ *it[0..2] , it[4].sort{ x,y -> sort_files(x,y) } ] }
	.set{ sampling }

/*|| add species and read length ||*/
/*||||*/
sampling
	.map{[
		*it[0..2],
		project.getSpecies( get_lims(it[1]) ),
		project.getRoughReadLength( get_lims(it[1]) ),
		it[3]
	]}
	.set{ sampling }

/*|| add rsem star index ||*/
/*||||*/
sampling
	.map{[
		*it[0..2],
		get_star_index( it[3] , it[4] ),
		it[5]
	]}
	.set{ sampling }

////////////////////
process sampled_star {
////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(genome),
			file(fastq) from sampling

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.Aligned.sortedByCoord.out.bam") into sampled_bam

	script:
		
		basename = flowcell + "_" + name
		cpus = task.cpus

		template "star.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------- SAMTOOLS -------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////////
process sampled_sort {
////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(bam) from sampled_bam

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.bam"),
			file("*.bai") \
			into \
				sampled_picard_rg,
				sampled_picard_rg_multiqc,
				sampled_sorted_bam,
				sampled_sorted_bam_multiqc

	script:

		basename = flowcell + "_" + name
		filename = basename + ".sorted.bam"

		template "samtools/sort_index.sh"
}

///////////////////////////////////////////////////////////////////////////////
// -------------------------------- PICARD --------------------------------- //
///////////////////////////////////////////////////////////////////////////////

//////////////////////////////
//process sampled_picard_group {
//////////////////////////////
//
//	tag { basename }
//
//	input:
//		set \
//			val(flowcell),
//			val(name),
//			val(single_end),
//			file(bam),
//			file(bai) from sampled_sorted_bam
//
//	output:
//		set \
//			val(flowcell),
//			val(name),
//			val(single_end),
//			file("*.rg.bam") \
//			into \
//				sampled_picard_rg,
//				sampled_picard_rg_multiqc
//
//	script:
//
//		basename = flowcell + "_" + name
//		filename = basename + ".rg.bam"
//		tmp_dirname = "tmp"
//
//		template "picard/group.sh"
//}

////////////////////////////////
process sampled_picard_duplicate {
////////////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(bam),
			file(bai) from sampled_picard_rg

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.dupmarked.bam") \
			into \
				sampled_picard_duplicate,
				sampled_picard_duplicate_multiqc

	script:

		basename = flowcell + "_" + name
		filename = basename + ".dupmarked.bam"
		metrics_filename = basename + ".marked_duplicates"
		tmp_dirname = "tmp"

		template "picard/duplicates.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------- SAMTOOLS -------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////////////////
process sampled_picard_index {
////////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(bam) from sampled_picard_duplicate

	output:
		// !! this output triggers a WARNING because we try to publish the bam
		// file a second time (the first time was with process
		// "sampled_picard_duplicate"
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(bam),
			file("*.bai") \
			into\
				sampled_bam_infer_experiment,
				sampled_bam_rnaseqc,
				sampled_picard_index_multiqc

	script:

		basename = flowcell + "_" + name

		template "samtools/index.sh"
}

///////////////////////////////////////////////////////////////////////////////
// --------------------------------- RSEQC --------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|| add bed file via species ||*/
/*||||||||||||||||||||||||*/
sampled_bam_infer_experiment
	.map{[
		*it[0..2],
		get_file( project.getSpecies( get_lims(it[1]) ) , "bed" ),
		*it[3..4]
	]}
	.set{ sampled_bam_infer_experiment }

////////////////////////////////
process sampled_infer_experiment {
////////////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(bed),
			file(bam),
			file(bai) from sampled_bam_infer_experiment
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.infer_experiment*") \
			into \
				sampled_infer_experiment,
				sampled_infer_experiment_multiqc
	
	script:

		basename = flowcell + "_" + name
		metrics_filename = basename + ".infer_experiment"

		template "rseqc/infer_experiment.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------- RNASEQC --------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|||||||||||||||*/
sampled_bam_rnaseqc
	.map{[
		*it[0..2],
		get_file(
			project.getSpecies( get_lims(it[1]) ),
			"fasta"
			),
		get_file(
			project.getSpecies( get_lims(it[1]) ),
			"rrna_list"
			),
		get_file(
			project.getSpecies( get_lims(it[1]) ),
			"rnaseqc_gtf"
			),
		*it[3..4]
	]}
	.set{ sampled_bam_rnaseqc }

///////////////////////
process sampled_rnaseqc {
///////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(fasta),
			val(rrna),
			val(gtf),
			file(bam),
			file(bai) from sampled_bam_rnaseqc

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("*.rnaseqc") into sampled_rnaseqc, sampled_rnaseqc_multiqc

	script:
		
		basename = flowcell + "_" + name
		metrics_filename = basename + ".rnaseqc"

		if (single_end) {

			single_end_arg = "-singleEnd"

		} else {

			single_end_arg = ""

		}

		template "rnaseqc.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------- MULTIQC --------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||||||||||||||||||||||*/
sampled_infer_experiment_multiqc
	.concat(sampled_rnaseqc_multiqc)
	.groupTuple(by:[0,1,2])
	.set{ to_sampled_multiqc }

///////////////////////
process sampled_multiqc {
///////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(metrics) from to_sampled_multiqc

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file("multiqc_data/multiqc_rseqc_infer_experiment.json"),
			file("multiqc_data/multiqc_rna_seqc.json") into sampled_multiqc

	script:

		basename = flowcell + "_" + name

		multiqc_data_format = "json"
		multiqc_template = "custom"
		multiqc_conf = multiqc_conf_file

		template "multiqc/sampled.sh"
}

///////////////////////////////////////////////////////////////////////////////
// --------------------------- STRANDEDNESS -------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////////
process strandedness {
////////////////////

	tag { basename }

	echo true
	executor "local"

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			file(rseqc_json),
			file(rnaseqc_json) from sampled_multiqc
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(strandedness) \
			into \
				strandedness,
				strandedness_align
	
	script:

		basename = flowcell + "_" + name
		prefix = single_end ? "se" : "pe"

		// get the json
		slurper = new JsonSlurper()
		rseqc = slurper.parse( new File( rseqc_json.toRealPath().toString() ) )
		metrics = rseqc[ basename + ".infer_experiment" ]

		// compute the metrics
		sense = metrics[ prefix + "_sense" ]
		antisense = metrics[ prefix + "_antisense" ]
		failed = metrics["failed"]
		normalised_sense = sense / (sense + antisense) * 100.0
		normalised_antisense = antisense / (sense + antisense) * 100.0

		// threshold
		if ( normalised_sense < 30 ) {
			strandedness = "reverse"
		} else if ( normalised_antisense > 80 ) {
			strandedness = "forward"
		} else {
			strandedness = "none"
		}

		"""
		"""
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ----------------------- ALIGNMENT SEQUENCING TYPE  ---------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||||||||||||||||||*/
rnaseq_cutadapt_align =
	Channel.create()

/*||||||||||||||||||||||||||||*/
non_rnaseq_cutadapt_align =
	Channel.create()

/*||||||||||*/
cutadapt_align
	.choice(rnaseq_cutadapt_align, non_rnaseq_cutadapt_align) {
		project.getSequencingType( get_lims(it[1]) ) == "rnaseq" ? 0 : 1
	}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// -------------------------- RNASEQ ALIGNMENT ----------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|||||||||||||||||*/
rnaseq_cutadapt_align
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": project.getSequencingType( get_lims(it[1]) ),
		"species": project.getSpecies( get_lims(it[1]) ),
		"read_length": project.getRoughReadLength( get_lims(it[1]) ),
		"fastq": it[3]
	]}
	.map{[
		it["flowcell"] + "_" + it["basename"],
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		get_rsem_star_index( it["species"] , it["read_length"] ),
		it["fastq"]
	]}
	.set{ cutadapt_align_strandedness }

/*||||||||||||||*/
strandedness_align
	.map{ [ it[0] + "_" + it[1] , it[3] ] }
	.combine(cutadapt_align_strandedness)
	.filter{ it[0] == it[2] }
	.map{ [ it[1] , *it[3..it.size()-1] ] }
	.set{ to_align_rsem_star }

/////////////////
process rsem_star {
/////////////////

	tag { basename }

	input:
		set \
			val(strandedness),
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(index),
			file(fastq) from to_align_rsem_star

	output:
		set \
			val(strandedness),
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.transcript.bam") into rsem_transcript
		set \
			val(strandedness),
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.results") into rsem_results
		set \
			val(strandedness),
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.stat") into rsem_stat
		set \
			val(strandedness),
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.STAR.genome.bam") \
			into \
				rsem_star_genome_sort,
				rsem_star_genome_multiqc

	script:

		cpus = task.cpus
		basename = flowcell + "_" + name

		if (single_end) {

			template "rsem/star/single_end.sh"

		} else {

			template "rsem/star/paired_end.sh"
		}
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------ NON RNASEQ ALIGNMENT --------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|||||||||||||||||||||*/
non_rnaseq_cutadapt_align
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": project.getSequencingType( get_lims(it[1]) ),
		"species": project.getSpecies( get_lims(it[1]) ),
		"read_length": project.getRoughReadLength( get_lims(it[1]) ),
		"fastq": it[3]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		get_star_index( it["species"] , it["read_length"] ),
		it["fastq"]
	]}
	.set{ to_align_star }

////////////
process star {
////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(genome),
			file(fastq) from to_align_star

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.out") into star_out
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.tab") into star_tab
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.Aligned.sortedByCoord.out.bam") into star_bam

	script:

		cpus = task.cpus
		basename = flowcell + "_" + name

		template "star.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------ SORTING AND INDEXING --------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|||||||||||||||||*/
rsem_star_genome_sort
	.map{ [ *it[1..it.size()-1] ] }
	.concat(star_bam)
	.set{ bam_to_sort }

///////////////////////////
process samtools_sort_index {
///////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam) from bam_to_sort

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.bam"),
			file("*.bai") \
			into \
				sorted_bam_picard,
				sorted_bam_multiqc

	script:

		basename = flowcell + "_" + name
		filename = basename + ".sorted.bam"

		template "samtools/sort_index.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------- PICARD ---------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// --------------------------- PICARD GROUP -------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////////
process picard_group {
////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam), file(bai) from sorted_bam_picard

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.bam") into picard_rg_duplicate, picard_rg_multiqc

	script:

		basename = flowcell + "_" + name

		tmp_dirname = "tmp"
		filename = basename + ".rg.bam"

		template "picard/group.sh"
}

///////////////////////////////////////////////////////////////////////////////
// -------------------------- PICARD DUPLICATE ----------------------------- //
///////////////////////////////////////////////////////////////////////////////

/////////////////////////
process picard_duplicates {
/////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam) from picard_rg_duplicate

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.marked_duplicates") into picard_duplicate
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.bam") \
			into \
				picard_bam_index,
				picard_bam_complexity,
				picard_bam_rnaseqmetrics,
				picard_bam_multimetrics,
				picard_bam_rnaseqc,
				picard_bam_qc_genes,
				picard_bam_transcripts

	script:

		basename = flowcell + "_" + name

		tmp_dirname = "tmp"
		filename = basename + ".dupmarked.bam"
		metrics_filename = basename + ".marked_duplicates"

		template "picard/duplicates.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ---------------------------- PICARD INDEX ------------------------------- //
///////////////////////////////////////////////////////////////////////////////

////////////////////
process picard_index {
////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam) from picard_bam_index

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam),
			file("*.bai") \
			into \
				picard_bai_rseqc,
				picard_bai_rnaseqc,
				picard_bai_tss,
				picard_bai_stats,
				picard_bai_idxstats,
				picard_bai_multiqc

	script:

		basename = flowcell + "_" + name

		template "samtools/index.sh"
}

///////////////////////////////////////////////////////////////////////////////
// --------------------------- PICARD COMPLEXITY --------------------------- //
///////////////////////////////////////////////////////////////////////////////

/////////////////////////
process picard_complexity {
/////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam) from picard_bam_complexity

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.complexity") into picard_complexity

	script:

		basename = flowcell + "_" + name

		tmp_dirname = "tmp"
		metrics_filename = basename + ".complexity"

		template "picard/complexity.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------- PICARD RNASEQMETRICS -------------------------- //
///////////////////////////////////////////////////////////////////////////////

//picard_bam_rnaseqmetrics
//	.map{[
//		"flowcell": it[0],
//		"basename": it[1],
//		"single_end": it[2],
//		"type": it[3],
//		"lims_id": get_lims(it[1]),
//		"species": project.getSpecies( get_lims(it[1]) ),
//		"bam": it[4]
//	]}
//	.println()

/*||||||||||||||||||||*/
picard_bam_rnaseqmetrics
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": it[3],
		"species": project.getSpecies( get_lims(it[1]) ),
		"bam": it[4]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		get_file( it["species"] , "refflat" ),
		get_file( it["species"] , "rrna_interval_list" ),
		it["bam"]
	]}
	.set{ to_picard_rnaseqmetrics }

////////////////////////////
process picard_rnaseqmetrics {
////////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(refflat),
			val(rrna_interval_list),
			file(bam) from to_picard_rnaseqmetrics
	
	when:
		type == "rnaseq"

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.rnaseqmetrics") into picard_rnaseqmetrics

	script:

		basename = flowcell + "_" + name

		tmp_dirname = "tmp"
		metrics_filename = basename + ".rnaseqmetrics"

		template "picard/rnaseqmetrics.sh"
}

///////////////////////////////////////////////////////////////////////////////
// -------------------------- RNASEQ MULTIMETRICS -------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|||||||||||||||||||*/
picard_bam_multimetrics
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": it[3],
		"species": project.getSpecies( get_lims(it[1]) ),
		"bam": it[4]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		get_file( it["species"] , "fasta" ),
		it["bam"]
	]}
	.set{ to_picard_multimetrics }

///////////////////////////
process picard_multimetrics {
///////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(fasta),
			file(bam) from to_picard_multimetrics

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.pdf") into picard_multimetrics_pdf
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*_metrics") into picard_multimetrics_metrics

	script:

		basename = flowcell + "_" + name

		tmp_dirname = "tmp"
		metrics_filename = basename + ".multimetrics"

		template "picard/multimetrics.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------ SAMTOOLS --------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

//////////////////////
process samtools_stats {
//////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam),
			file(bai) from picard_bai_stats

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.stats") into stats

	script:

		basename = flowcell + "_" + name
		filename = basename + ".stats"

		template "samtools/stats.sh"
}

/////////////////////////
process samtools_idxstats {
/////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(bam),
			file(bai) from picard_bai_idxstats

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.idxstats") \
			into \
				idxstats_chrom,
				idxstats_multiqc

	script:

		basename = flowcell + "_" + name
		filename = basename + ".idxstats"

		template "samtools/idxstats.sh"
}

//////////////////
process chromosome {
//////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file(idxstats) from idxstats_chrom

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*_chrom.json") into chrom
	
	shell:

		basename = flowcell + "_" + name
		filename = basename + "_chrom.json"

		"""
		python ${py_dir}/reads_per_chrom.py ${idxstats} > ${filename}
		"""
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// -------------------------------- RSEQC ---------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||||||*/
picard_bai_rseqc
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": it[3],
		"rough_read_length": project.getRoughReadLength( get_lims(it[1]) ),
		"bed": get_file( project.getSpecies( get_lims(it[1]) ), "bed" ),
		"bam": it[4],
		"bai": it[5]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		it["rough_read_length"],
		it["bed"],
		it["bam"],
		it["bai"]
	]}
	.into{
		picard_bai_infer_experiment;
		picard_bai_junction_annotation;
		picard_bai_junction_saturation;
		picard_bai_mismatch_profile;
		picard_bai_read_distribution;
		picard_bai_transcript_integrity;
	}

////////////////////////
process infer_experiment {
////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(rough_read_length),
			val(bed),
			file(bam),
			file(bai) from picard_bai_infer_experiment
	
	when:
		type == "rnaseq"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.infer_experiment*") into infer_experiment
	
	script:

		basename = flowcell + "_" + name
		metrics_filename = basename + ".infer_experiment"

		template "rseqc/infer_experiment.sh"
}

///////////////////////////
process junction_annotation {
///////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(rough_read_length),
			val(bed),
			file(bam),
			file(bai) from picard_bai_junction_annotation

	when:
		type == "rnaseq"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.junction_annotation*") into junction_annotation
	
	script:

		basename = flowcell + "_" + name
		metrics_filename = basename + ".junction_annotation"

		template "rseqc/junction_annotation.sh"
}

///////////////////////////
process junction_saturation {
///////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(rough_read_length),
			val(bed),
			file(bam),
			file(bai) from picard_bai_junction_saturation

	when:
		type == "rnaseq"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.junction_saturation*") into junction_saturation

	script:

		basename = flowcell + "_" + name
		metrics_filename = basename + ".junction_saturation"

		template "rseqc/junction_saturation.sh"
}

////////////////////////
process mismatch_profile {
////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(rough_read_length),
			val(bed),
			file(bam),
			file(bai) from picard_bai_mismatch_profile

	when:
		type == "rnaseq"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.mismatch_profile*") into mismatch_profile

	script:

		basename = flowcell + "_" + name
		metrics_filename = basename + ".mismatch_profile"

		template "rseqc/mismatch_profile.sh"
}

/////////////////////////
process read_distribution {
/////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(rough_read_length),
			val(bed),
			file(bam),
			file(bai) from picard_bai_read_distribution

	when:
		type == "rnaseq"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.read_distribution*") into read_distribution
	
	script:

		basename = flowcell + "_" + name
		metrics_filename = basename + ".read_distribution"

		template "rseqc/read_distribution.sh"
}

///////////////////////////////////
process transcript_integrity_number {
///////////////////////////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(rough_read_length),
			val(bed),
			file(bam),
			file(bai) from picard_bai_transcript_integrity

	when:
		type == "rnaseq"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.{xls,txt}") into transcript_integrity

	script:

		basename = flowcell + "_" + name

		template "rseqc/transcript_integrity_number.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------ RNASEQC ---------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||||||||*/
picard_bai_rnaseqc
	.map{[
		*it[0..3],
		get_file(
			project.getSpecies( get_lims(it[1]) ),
			"fasta"
			),
		get_file(
			project.getSpecies( get_lims(it[1]) ),
			"rrna_list"
			),
		get_file(
			project.getSpecies( get_lims(it[1]) ),
			"rnaseqc_gtf"
			),
		*it[4..5]
	]}
	.set{ to_rnaseqc }

///////////////
process rnaseqc {
///////////////

	tag { basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(fasta),
			val(rrna),
			val(gtf),
			file(bam),
			file(bai) from to_rnaseqc

	when:
		type == "rnaseq"

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*.rnaseqc") into rnaseqc

	script:
		
		basename = flowcell + "_" + name
		metrics_filename = basename + ".rnaseqc"

		if (single_end) {

			single_end_arg = "-singleEnd"

		} else {

			single_end_arg = ""

		}

		template "rnaseqc.sh"
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ---------------------------- TSS COVERAGE ------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||||*/
picard_bai_tss
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": it[3],
		"gtf": get_file( project.getSpecies( get_lims(it[1]) ), "gtf" ),
		"bam": it[4],
		"bai": it[5]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		it["gtf"],
		it["bam"],
		it["bai"]
	]}
	.set{ to_tss }

///////////
process tss {
///////////

	tag{ basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(gtf),
			file(bam),
			file(bai) from to_tss
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*_tss.json") into tss
	
	shell:

		basename = flowcell + "_" + name

		"""
		python ${py_dir}/tss_coverage.py \
			--bam ${bam} \
			--gtf ${gtf} \
			--win 3000 \
			--size 200 \
			--name ${basename} \
			> ${basename}_tss.json
		"""
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ----------------------- QC GENES STRANDEDNESS --------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*|||||||||||||||*/
picard_bam_qc_genes
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": it[3],
		"species": project.getSpecies( get_lims(it[1]) ),
		"bam": it[4]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		it["species"],
		get_file( it["species"] , "bed" ),
		Paths.get(qc_genes_dir, "human_rrna_ids.txt"),
		Paths.get(qc_genes_dir, "human_mtrrna_ids.txt"),
		Paths.get(qc_genes_dir, "human_globin_ids.txt"),
		it["bam"]
	]}
	.set{ to_qc_genes }

////////////////
process qc_genes {
////////////////

	tag{ basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(species),
			val(bed),
			val(rrna),
			val(mtrrna),
			val(globin),
			file(bam) from to_qc_genes
	
	when:
		species.toLowerCase() == "homo sapiens"
	
	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*_qc_genes.json") into qc_genes
	
	shell:

		basename = flowcell + "_" + name

		"""
		${cpp_dir}/qc_directionality \
			${bam} \
			${bed} \
			${rrna} \
			${mtrrna} \
			${globin} \
			> ${basename}_qc_genes.json
		"""
}

/*||||||||||||||||||*/
picard_bam_transcripts
	.map{[
		"flowcell": it[0],
		"basename": it[1],
		"single_end": it[2],
		"type": it[3],
		"bed": get_file( project.getSpecies( get_lims(it[1]) ), "bed" ),
		"bam": it[4]
	]}
	.map{[
		it["flowcell"],
		it["basename"],
		it["single_end"],
		it["type"],
		it["bed"],
		it["bam"]
	]}
	.set{ to_transcripts }

///////////////////
process transcripts {
///////////////////

	tag{ basename }

	input:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			val(bed),
			file(bam) from to_transcripts
	
	when:
		type == "rnaseq"

	output:
		set \
			val(flowcell),
			val(name),
			val(single_end),
			val(type),
			file("*_transcripts.json") into transcripts
	
	shell:

		basename = flowcell + "_" + name

		"""
		${cpp_dir}/transcript_directionality \
			${bam} \
			${bed} \
			> ${basename}_transcripts.json
		"""
}

///////////////////////////////////////////////////////////////////////////////
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------- MULTIQC --------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
///////////////////////////////////////////////////////////////////////////////

/*||||||||*/
info_multiqc =
	Channel
		.from( project.getFastq() )
		.map{ get_experiment_info(it) }
		.flatten()
		.collate(2)
		.groupTuple()
		.map{ [ it[0] , it[1].toUnique().sort().join(",") ] }
		.toList()
		.map{ it.collectEntries{ x -> [ (x[0]): x[1] ] } }
		.collect()
		.get()
		.getAt(0)

////////////////////
process multiqc_conf {
////////////////////

	output:
		file "multiqc_conf.yml" into multiqc_run_conf

	shell:

		def filename = "multiqc_conf.yml"

		def logo = Paths.get(
						workflow.projectDir.toString(),
						"conf",
						"logo.jpg"
						).toString()

		// additional details
		info_multiqc["Contact e-mail"] = "bioinformatics@crick.ac.uk"

		def conf = [:]
		conf["custom_logo"] = logo
		conf["report_header_info"] = info_multiqc

		def json = JsonOutput.toJson(conf)

		def header = []
		info_multiqc.each{ k, v ->
			s = "   - " + k + ": " + v
			header.add(s)
		}
		def header_string = header.join("\\n")

		"""
		cat ${multiqc_conf_file} > ${filename}
		echo "custom_logo: ${logo}" >> ${filename}
		echo "report_header_info:" >> ${filename}
		echo -e "${header_string}" >> ${filename}
		"""
}

////////////
process info {
////////////

	output:
		file "*.json" into experiment_multiqc_info

	shell:

		def info = [
			protocol: "test",
			strandedness: "test",
			read_length: "test",
			species: "test"
		]

		def json = JsonOutput.toJson(info)

		"""
		echo -E '${json}' > experiment_info_transcripts.json
		echo -E '${json}' > experiment_info_qc_genes.json
		"""
}

///////////////
process multiqc {
///////////////

	input:
		
		file multiqc_conf from multiqc_run_conf
		file experiment_multiqc_conf from experiment_multiqc_info
		
		file("*") from cutadapt_log.map{x->x[3]}.collect().ifEmpty([])
		file("*") from cutadapt_multiqc.map{x->x[3]}.collect().ifEmpty([])

		file("*") from fastqc.map{x->x[2]}.collect().ifEmpty([])
		file("*") from fastqc_cutadapt.map{x->x[4]}.collect().ifEmpty([])

		file("*") from fscreen_html.map{x->x[3]}.collect().ifEmpty([])
		file("*") from fscreen_txt.map{x->x[3]}.collect().ifEmpty([])

		file("*") from rsem_transcript.map{x->x[5]}.collect().ifEmpty([])
		file("*") from rsem_results.map{x->x[5]}.collect().ifEmpty([])
		file("*") from rsem_stat.map{x->x[5]}.collect().ifEmpty([])
		file("*") from \
			rsem_star_genome_multiqc.map{x->x[5]}.collect().ifEmpty([])

		file("*") from picard_rg_multiqc.map{x->x[4]}.collect().ifEmpty([])
		file("*") from picard_duplicate.map{x->x[4]}.collect().ifEmpty([])
		file("*") from picard_complexity.map{x->x[4]}.collect().ifEmpty([])
		file("*") from picard_rnaseqmetrics.map{x->x[4]}.collect().ifEmpty([])
		file("*") from \
			picard_multimetrics_metrics.map{x->x[4]}.collect().ifEmpty([])

		file("*") from stats.map{x->x[4]}.collect().ifEmpty([])
		file("*") from idxstats_multiqc.map{x->x[4]}.collect().ifEmpty([])

		file("*") from infer_experiment.map{x->x[4]}.collect().ifEmpty([])
		file("*") from junction_annotation.map{x->x[4]}.collect().ifEmpty([])
		file("*") from junction_saturation.map{x->x[4]}.collect().ifEmpty([])
		file("*") from mismatch_profile.map{x->x[4]}.collect().ifEmpty([])
		file("*") from read_distribution.map{x->x[4]}.collect().ifEmpty([])
		file("*") from transcript_integrity.map{x->x[4]}.collect().ifEmpty([])

		file("*") from rnaseqc.map{x->x[4]}.collect().ifEmpty([])

		file("*") from chrom.map{x->x[4]}.collect().ifEmpty([])
		file("*") from tss.map{x->x[4]}.collect().ifEmpty([])
		file("*") from qc_genes.map{x->x[4]}.collect().ifEmpty([])
		file("*") from transcripts.map{x->x[4]}.collect().ifEmpty([])

	output:
		file "multiqc_data" into multiqc_data
		file "multiqc_report.html" into multiqc_report

	script:

		multiqc_data_format = "json"
		multiqc_template = "custom"

		template "multiqc/general.sh"
}

