
// c header
#include <stdlib.h>
#include <limits.h>
#include <zlib.h>

// cpp header
#include <iostream>
#include <string>
#include <vector>
#include <set>

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
	std::cout << "connard !" << std::endl;
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
typedef seqan::IntervalAndCargo<int, seqan::CharString> TInterval;
typedef seqan::String<TInterval> StringTInterval;
typedef std::map< std::string, StringTInterval > GRange;
typedef seqan::IntervalTree<int, seqan::CharString> Tree;
typedef std::map< std::string, Tree> GTree;
typedef std::vector<std::string> Strands;
typedef std::map< std::string, Strands> TStrands;
typedef std::map< std::string, std::map< std::string, double >> TCount;

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

	// the path of the bam file and the path of the bed file
	std::string bam_path = abspath(argv[1]);
	std::string bed_path = abspath(argv[2]);

	// for file reading
	seqan::BamFileIn bam;
	seqan::BamHeader hdr;
	seqan::BamAlignmentRecord bam_rec;
	seqan::BedFileIn bed;
	seqan::BedRecord<seqan::Bed6> bed_rec;

	// strandedness strings
	std::vector<std::string> paired_sense_strings{"1++", "1--", "2+-", "2-+"};
	std::vector<std::string>
		paired_antisense_strings{"1+-", "1-+", "2++", "2--"};
	std::vector<std::string> single_sense_strings{"++", "--"};
	std::vector<std::string> single_antisense_strings{"+-", "-+"};

	////////////////////////////////////////////////////////////////////////////
	// CREATE A LIST OF THE TRANSCRIPT WITH THEIR STRAND
	
	// open the bed file
	if ( ! open(bed, bed_path.c_str()) )
	{
		std::cerr << "Error: cannot open " << bed_path << std::endl;
	}

	// strand of reads for a transcript
	TStrands mapping;

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

		// initialise the mapping
		Strands s;
		mapping[name] = s;
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

	// the required sam flags
	int flag_paired = 1;  // 0x1
	int flag_unmap = 4; // 0x4
	int flag_rev = 16; // 0x10
	int flag_read1 = 64; // 0x40
	int flag_read2 = 128; // 0x40
	int flag_secondary = 256; // 0x100
	int flag_qcfail = 512; // 0x200
	int flag_dup = 1024; // 0x400

	// here we go!
	while ( ! seqan::atEnd(bam) )
	{

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

		// the strand for each transcript overlap
		std::vector<std::string> strands;
		std::vector<std::string> matched_transcripts;
		seqan::String<seqan::CharString> results;
		seqan::findIntervals(results, trees[ref], begin, end);
		for (unsigned i=0; i<seqan::length(results); i++)
		{
			std::string transcript = seqan::toCString(results[i]);

			// strand
			std::string strand( 1, transcript.back() );
			strands.push_back(strand);

			// name
			std::string matched = transcript.substr(0, transcript.size()-1);
			matched_transcripts.push_back(matched);
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

			// strand determination
			std::string directionality;
			if (
					std::find(
						paired_sense_strings.begin(), 
						paired_sense_strings.end(),
						key) != paired_sense_strings.end()
				) {

					directionality = "+";

			} else if (
					std::find(
						paired_antisense_strings.begin(), 
						paired_antisense_strings.end(),
						key) != paired_antisense_strings.end()
					) {

					directionality = "-";

			} else {

					directionality = "*";
			}

			// mapping
			for (
					auto it = matched_transcripts.begin();
					it != matched_transcripts.end();
					++it
					) {
				mapping[*it].push_back(directionality);
			}

		} else {

			// key
			std::string key = reverse + strand;

			// strand determination
			std::string directionality;
			if (
					std::find(
						single_sense_strings.begin(), 
						single_sense_strings.end(),
						key) != single_sense_strings.end()
				) {

					directionality = "+";

			} else if (
					std::find(
						single_antisense_strings.begin(), 
						single_antisense_strings.end(),
						key) != single_antisense_strings.end()
					) {

					directionality = "-";

			} else {

					directionality = "*";
			}

			// mapping
			for (
					auto it = matched_transcripts.begin();
					it != matched_transcripts.end();
					++it
					) {
				mapping[*it].push_back(directionality);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////

	// counting for each transcript
	TCount transcript_count;
	for (auto it=mapping.begin(); it!=mapping.end(); ++it)
	{
		if (it->second.size() > 0) {
			int sense = std::count(it->second.begin(), it->second.end(), "+");
			int antisense = std::count(it->second.begin(), it->second.end(), "-");
			int undeter = std::count(it->second.begin(), it->second.end(), "*");
			std::map< std::string, double > count;
			count["sense"] = (double)sense / (double)it->second.size();
			count["antisense"] = (double)antisense / (double)it->second.size();
			count["undetermined"] = (double)undeter / (double)it->second.size();
			transcript_count[it->first] = count;
		}
	}

	// print
	std::cout << "{" << std::endl;
	for (auto i=transcript_count.begin(); i!=transcript_count.end();)
	{
		std::cout << "\t" << '"' << i->first << '"' << ":{";
		for (auto j=i->second.begin(); j!=i->second.end();)
		{
			std::cout << '"' << j->first << '"' << ":" << j->second;
			++j;

			if ( j != i->second.end() ) {
				std::cout << ",";
			}
		}
		std::cout << "}";
		++i;
		if ( i != transcript_count.end() ) {
			std::cout << ",";
		}
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;

	return 0;
}

