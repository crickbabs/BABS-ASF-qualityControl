
/*
 * nourdine.bah@crick.ac.uk
 *
 * this code implements the approaches of RSeQC and RNA-SeQC to quantify the
 * sense and antisense reads
 *
 * http://rseqc.sourceforge.net/#infer-experiment-py
 * http://archive.broadinstitute.org/cancer/cga/rna-seqc
 *
 * thanks to them!
 */

// c header
#include <stdlib.h>
#include <limits.h>
#include <zlib.h>

// cpp header
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>

// seqan header
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/bed_io.h>
#include <seqan/misc/interval_tree.h>

///////////////////////////////////////////////////////////////////////////////
// usage
void help_message()
{
	std::cout << "you idiot!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// absolute path
std::string abspath(char *path)
{
	char *full_path = realpath(path, NULL);
	std::string abs_path = full_path;
	free(full_path);
	return abs_path;
}

///////////////////////////////////////////////////////////////////////////////
// couting methods

/* ----------------- */
int get_strand_counting(
		std::vector<std::string> keys,
		std::map<std::string, int> counting
		)
{
	int sum = 0;
	std::vector<std::string>::iterator key;
	for (auto key=keys.begin(); key!=keys.end(); ++key)
	{
		if (counting.count(*key) != 0) {
			sum += counting[*key];
		}
	}
	return sum;
}

/* ----------------------------------------------------- */
int get_single_sense(std::map<std::string, int> counting)
{
	std::vector<std::string> keys = {"++", "--"};
	return get_strand_counting(keys, counting);
}

/* ----------------------------------------------------- */
int get_single_antisense(std::map<std::string, int> counting)
{
	std::vector<std::string> keys = {"+-", "-+"};
	return get_strand_counting(keys, counting);
}

/* ----------------------------------------------------- */
int get_paired_sense(std::map<std::string, int> counting)
{
	std::vector<std::string> keys = {"1++", "1--", "2+-", "2-+"};
	return get_strand_counting(keys, counting);
}

/* --------------------------------------------------------- */
int get_paired_antisense(std::map<std::string, int> counting)
{
	std::vector<std::string> keys = {"1+-", "1-+", "2++", "2--"};
	return get_strand_counting(keys, counting);
}

///////////////////////////////////////////////////////////////////////////////
typedef seqan::IntervalAndCargo<int, seqan::CharString> TInterval;
typedef seqan::String<TInterval> StringTInterval;
typedef std::map< std::string, StringTInterval > GRange;
typedef seqan::IntervalTree<int, seqan::CharString> Tree;
typedef std::map< std::string, Tree> GTree;
typedef std::map< std::string, double> Result;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	// check if zlib in the path
	#if ! SEQAN_HAS_ZLIB
		std::cout << "zlib not found" << std::endl;
		exit(1);
	#endif

	// the number of arguments
	if (argc<3) {
		help_message();
		exit(1);
	}

	// the path of the bam file and the path of the bed files
	std::string bam_path = abspath(argv[1]);
	std::string bed_path = abspath(argv[2]);
	std::string rrna_txt_path = abspath(argv[3]);
	std::string mtrrna_txt_path = abspath(argv[4]);
	std::string globin_txt_path = abspath(argv[5]);

	// the files contain a list of transcript ids
	std::ifstream rrna_txt(rrna_txt_path);
	std::ifstream mtrrna_txt(mtrrna_txt_path);
	std::ifstream globin_txt(globin_txt_path);

	// the transcript ids in a vector in order to be able to search
	std::vector<std::string> rrna_ids;
	std::vector<std::string> mtrrna_ids;
	std::vector<std::string> globin_ids;
	std::string line;
	while ( std::getline(rrna_txt, line) ) { rrna_ids.push_back(line); }
	while ( std::getline(mtrrna_txt, line) ) { mtrrna_ids.push_back(line); }
	while ( std::getline(globin_txt, line) ) { globin_ids.push_back(line); }

	// for file reading
	seqan::BamFileIn bam;
	seqan::BamHeader hdr;
	seqan::BamAlignmentRecord bam_rec;
	seqan::BedFileIn bed;
	seqan::BedRecord<seqan::Bed6> bed_rec;

	////////////////////////////////////////////////////////////////////////////
	// CREATE A LIST OF THE TRANSCRIPT WITH THEIR STRAND
	
	// open the bed file
	if ( ! open(bed, bed_path.c_str()) )
	{
		std::cerr << "Error: cannot open " << bed_path << std::endl;
	}

	// the genomic regions of the bed file
	GRange granges;

	while ( ! seqan::atEnd(bed) )
	{
		seqan::readRecord(bed_rec, bed);

		// transcript features
		std::string name = seqan::toCString(bed_rec.name);
		std::string ref = seqan::toCString(bed_rec.ref);
		int begin = bed_rec.beginPos;
		int end = bed_rec.endPos;
		std::string strand(1, bed_rec.strand);

		// an interval list for each chromosome
		if (granges.count(ref) == 0) {
			seqan::String<TInterval> intervals;
			granges[ref] = intervals;
		}

		// add the interval to the string
		seqan::appendValue( granges[ref], TInterval(begin, end, name+strand) );
	}

	// create a tree from the intervals for each chromosome
	GTree trees;
	GRange::iterator it;
	for (auto it=granges.begin(); it!=granges.end(); ++it)
	{
		Tree tree(it->second);
		trees[it->first] = tree;
	}

	////////////////////////////////////////////////////////////////////////////
	// GET STRAND FOR EACH ENTRY OF THE BAM FILE

	// open the bam file
	if ( ! open(bam, bam_path.c_str()) )
	{
		std::cerr << "Error: cannot open " << bam_path << std::endl;
		exit(1);
	}

	// get the context
	typedef seqan::FormattedFileContext<
		seqan::BamFileIn,
		void
		>::Type t_bam_context;
	t_bam_context const &bam_context = context(bam);

	// dunno why but seqan cannot read the entries if it doesn't read the header
	// first...
	seqan::readHeader(hdr, bam);

	// useful to check the progression
	int total = 0;
	int success = 0;
	int mapped = 0;
	//int total_offset = 1e6;

	// the required sam flags
	int flag_paired = 1;  // 0x1
	int flag_unmap = 4; // 0x4
	int flag_rev = 16; // 0x10
	int flag_read1 = 64; // 0x40
	int flag_read2 = 128; // 0x40
	int flag_secondary = 256; // 0x100
	int flag_qcfail = 512; // 0x200
	int flag_dup = 1024; // 0x400

	// counting
	std::map<std::string, int> single;
	std::map<std::string, int> paired;
	std::map<std::string, int> single_rrna;
	std::map<std::string, int> paired_rrna;
	std::map<std::string, int> single_mtrrna;
	std::map<std::string, int> paired_mtrrna;
	std::map<std::string, int> single_globin;
	std::map<std::string, int> paired_globin;
	std::map<std::string, int> single_nonqc;
	std::map<std::string, int> paired_nonqc;

	// here we go!
	while ( ! seqan::atEnd(bam) )
	{

		total++;
		/*
		if (total % (int)1e6 == 0) {
			std::cout << total << std::endl;
		}
		// */

		// the current read
		seqan::readRecord(bam_rec, bam);
		int flag = bam_rec.flag;

		// only the uniquely mapped reads
		if (
				(flag & flag_unmap) ||
				(flag & flag_secondary) ||
				(flag & flag_qcfail) ||
				(flag & flag_dup)
				)
		{ continue; }

		// skip multi mapped reads
		bool multi_mapped = false;
		seqan::BamTagsDict tags(bam_rec.tags);
		for (unsigned id=0; id<seqan::length(tags); id++)
		{
			if (seqan::getTagKey(tags, id) == "NH") {

				int nh = 0;

				if ( !seqan::extractTagValue(nh, tags, id) ) {
					std::cout << "Problem: cannot extract an NH tag value."
						<< std::endl;
					exit(1);
				}

				if ( nh > 1 ) {
					multi_mapped = true;
				}
			}
		}
		if (multi_mapped) {
			continue;
		} else {
			mapped++;
		}

		/////////////////////////////////////////////////////////////////////////
		// read features

		int rid = bam_rec.rID;
		seqan::String<char> cname = seqan::contigNames(bam_context)[rid];
		std::string ref = seqan::toCString(cname);

		seqan::String<char> cseq = bam_rec.seq;
		std::string sequence = seqan::toCString(cseq);
		int length = sequence.length();

		int begin = bam_rec.beginPos;
		int end = begin + length;

		std::string reverse = flag & flag_rev ? "-" : "+";
		/////////////////////////////////////////////////////////////////////////

		// qc genes flags
		bool is_rrna = false;
		bool is_mtrrna = false;
		bool is_globin = false;

		// the strand for each transcript overlap
		std::vector<std::string> strands;
		seqan::String<seqan::CharString> results;
		seqan::findIntervals(results, trees[ref], begin, end);
		for (unsigned i=0; i<seqan::length(results); i++)
		{
			// we need to check if of all the strands are the same
			std::string transcript = seqan::toCString(results[i]);
			std::string strand( 1, transcript.back() );
			strands.push_back(strand); 

			// we need to check if the read is overlapping a qc gene
			std::string transcript_id = transcript.substr(0, transcript.size()-1);
			if ( std::find(rrna_ids.begin(), rrna_ids.end(), transcript_id)
					!=
					rrna_ids.end() ) {
				is_rrna = true;
			}
			if ( std::find(mtrrna_ids.begin(), mtrrna_ids.end(), transcript_id)
					!=
					mtrrna_ids.end() ) {
				is_mtrrna = true;
			}
			if ( std::find(globin_ids.begin(), globin_ids.end(), transcript_id)
					!=
					globin_ids.end() ) {
				is_globin = true;
			}
		}

		// all overlaps have to have the same strand
		std::string strand;
		std::set<std::string> uniq_strands(strands.begin(), strands.end());
		if (uniq_strands.size() == 1) {
			std::vector<std::string> ss;
			ss.assign(uniq_strands.begin(), uniq_strands.end());
			strand = ss[0];
		} else {
			continue;
		}

		/////////////////////////////////////////////////////////////////////////
		// same approach as in infer_experiment.py of RSeQC
		if (flag & flag_paired) {

			// key
			std::string read;
			if (flag & flag_read1) { read = "1"; }
			if (flag & flag_read2) { read = "2"; }
			std::string key = read + reverse + strand;

			// set keys to 0
			if (paired.count(key) == 0) { paired[key] = 0; }
			if (paired_rrna.count(key) == 0) { paired_rrna[key] = 0; }
			if (paired_mtrrna.count(key) == 0) { paired_mtrrna[key] = 0; }
			if (paired_globin.count(key) == 0) { paired_globin[key] = 0; }
			if (paired_nonqc.count(key) == 0) { paired_nonqc[key] = 0; }

			// counting
			success++;
			paired[key]++;

			// qc
			if (is_rrna) {
				paired_rrna[key]++;
			}
			else if (is_mtrrna) {
				paired_mtrrna[key]++;
			}
			else if (is_globin) {
				paired_globin[key]++;
			} else {
				paired_nonqc[key]++;
			}

		} else {

			// key
			std::string key = reverse + strand;

			// set keys to 0
			if (single.count(key) == 0) { single[key] = 0; }
			if (single_rrna.count(key) == 0) { single_rrna[key] = 0; }
			if (single_mtrrna.count(key) == 0) { single_mtrrna[key] = 0; }
			if (single_globin.count(key) == 0) { single_globin[key] = 0; }
			if (single_nonqc.count(key) == 0) { single_nonqc[key] = 0; }

			// counting
			success++;
			single[key]++;

			// qc
			if (is_rrna) {
				single_rrna[key]++;
			}
			if (is_mtrrna) {
				single_mtrrna[key]++;
			}
			if (is_globin) {
				single_globin[key]++;
			} else {
				single_nonqc[key]++;
			}
		}
		/////////////////////////////////////////////////////////////////////////

		// todo: maybe an alternative : need to check !
		//bool reverse = flag & flag_rev;
		//if ( reverse && sense ) {
		//	positive++;
		//} else if ( !reverse && !sense ) {
		//	positive++;
		//} else if ( reverse && !sense ) {
		//	negative++;
		//} else if ( !reverse && sense ) {
		//	negative++;
		//}
	}

	//std::map<std::string, int>::iterator iit;
	//for (auto iit=paired.begin(); iit!=paired.end(); ++iit) {
	//	std::cout << iit->first << ":" << iit->second << std::endl;
	//}
	//for (auto iit=single.begin(); iit!=single.end(); ++iit) {
	//	std::cout << iit->first << ":" << iit->second << std::endl;
	//}
	
	////////////////////////////////////////////////////////////////////////////
	// EXPORT RESULTS


	std::string protocol;

	Result result;
	result["total"] = total;
	result["mapped"] = mapped;
	result["success"] = success;


	if ( paired.size()>0 && single.size()==0 ) {

		protocol = "paired end";

		// the sum
		int paired_success = 0;
		std::map<std::string, int>::iterator it;
		for (auto it=paired.begin(); it!=paired.end(); ++it)
		{
			paired_success += it->second;
		}

		// the choice of normalisation
		//int norm_count = total;
		//int norm_count = success;
		//int norm_count = paired_success;
		int norm_count = mapped;

		// counting
		int sense = get_paired_sense(paired);
		int antisense = get_paired_antisense(paired);
		int sense_rrna = get_paired_sense(paired_rrna);
		int antisense_rrna = get_paired_antisense(paired_rrna);
		int sense_mtrrna = get_paired_sense(paired_mtrrna);
		int antisense_mtrrna = get_paired_antisense(paired_mtrrna);
		int sense_globin = get_paired_sense(paired_globin);
		int antisense_globin = get_paired_antisense(paired_globin);
		int sense_nonqc = get_paired_sense(paired_nonqc);
		int antisense_nonqc = get_paired_antisense(paired_nonqc);

		// percentages
		double sense_percent = (double)sense / (double)norm_count;
		double antisense_percent = (double)antisense / (double)norm_count;
		double failed = 1 - sense_percent - antisense_percent;

		result["paired success"] = (double)paired_success;

		result["sense percent"] = sense_percent;
		result["antisense percent"] = antisense_percent;
		result["failed percent"] = failed;

		result["sense rRNA"] = (double)sense_rrna;
		result["antisense rRNA"] = (double)antisense_rrna;
		result["sense mt-rRNA"] = (double)sense_mtrrna;
		result["antisense mt-rRNA"] = (double)antisense_mtrrna;
		result["sense globin"] = (double)sense_globin;
		result["antisense globin"] = (double)antisense_globin;
		result["sense non QC"] = (double)sense_nonqc;
		result["antisense non QC"] = (double)antisense_nonqc;


	} else if ( paired.size()==0 && single.size()>0 ) {

		protocol = "single end";

		// the sum
		int single_success = 0;
		std::map<std::string, int>::iterator it;
		for (auto it=single.begin(); it!=single.end(); ++it)
		{
			single_success += it->second;
		}

		// the choice of normalisation
		//int norm_count = total;
		//int norm_count = success;
		//int norm_count = single_success;
		int norm_count = mapped;

		// counting
		int sense = get_single_sense(single);
		int antisense = get_single_antisense(single);
		int sense_rrna = get_single_sense(single_rrna);
		int antisense_rrna = get_single_antisense(single_rrna);
		int sense_mtrrna = get_single_sense(single_mtrrna);
		int antisense_mtrrna = get_single_antisense(single_mtrrna);
		int sense_globin = get_single_sense(single_globin);
		int antisense_globin = get_single_antisense(single_globin);
		int sense_nonqc = get_single_sense(single_nonqc);
		int antisense_nonqc = get_single_antisense(single_nonqc);
		
		// percentage
		double sense_percent = (double)sense / (double)norm_count;
		double antisense_percent = (double)antisense / (double)norm_count;
		double failed = 1 - sense_percent - antisense_percent;

		result["single success"] = (double)single_success;

		result["sense percent"] = sense_percent;
		result["antisense percent"] = antisense_percent;
		result["failed percent"] = failed;

		result["sense rRNA"] = (double)sense_rrna;
		result["antisense rRNA"] = (double)antisense_rrna;
		result["sense mt-rRNA"] = (double)sense_mtrrna;
		result["antisense mt-rRNA"] = (double)antisense_mtrrna;
		result["sense globin"] = (double)sense_globin;
		result["antisense globin"] = (double)antisense_globin;
		result["sense non QC"] = (double)sense_nonqc;
		result["antisense non QC"] = (double)antisense_nonqc;

	} else {

		protocol = "mixed";

	}

	// json output
	std::cout << "{" << std::endl;
	std::cout
		<< "\t" << '"' << "protocol" << '"' << ':'
		<< '"' << protocol << '"' << ','
		<< std::endl;
	Result::iterator res;
	for (auto res=result.begin(); res!=result.end();)
	{
		std::cout
			<< "\t" << '"'
			<< res->first
			<< '"' << ':'
			<< res->second;

		++res;

		if ( res != result.end() ) {
			std::cout << ',';
		}

		std::cout << std::endl;
	}

	std::cout << "}" << std::endl;

	return 0;
}

