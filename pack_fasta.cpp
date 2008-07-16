#ifdef PACK_FASTA_MAIN
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include "packed_io.h"
#include "tokenize.h"

using namespace std;
using namespace seqan;

static bool dna     = true;
static bool rna     = false;
static bool pro     = false;
static bool verbose = false;
static bool pack    = true;  // pack or unpack?

static void print_usage() {
	cout << "Usage: PackFasta [-d|-r|-p] <in_fasta_file> <out_packed_file>" << endl;
	cout << "    -d     input file uses 4-char DNA character set" << endl;
	cout << "    -r     input file uses 4-char RNA character set" << endl;
	cout << "    -p     input file uses 20+-char protein character set" << endl;
	cout << "    -u     unpack (default: pack)" << endl;
	cout << "    -v     verbose output" << endl;
}

int main(int argc, char **argv) {
	string infile;
	vector<string> infiles;
	string outfile;
	const char *short_options = "uvdrp";
	int next_option;
	do { 
		next_option = getopt(argc, argv, short_options);
		switch (next_option) {
			case 'd': /* DNA */
				dna = true;
				rna = false;
				pro = false;
				break;
			case 'r': /* RNA */
				dna = false;
				rna = true;
				pro = false;
	   		case 'p': /* Protein */
				dna = false;
				rna = false;
				pro = true;
				break;
	   		case 'v': /* verbose */
				verbose = true;
				break;
	   		case 'u': /* unpack */
				pack = false;
				break;
			case -1: /* Done with options. */
				break;
			default: 
				print_usage();
				return 1;
		}
	} while(next_option != -1);

	if(optind >= argc) {
		print_usage();
		return 1;
	}
	infile = argv[optind++];

	tokenize(infile, ",", infiles);
	if(infiles.size() < 1) {
		cerr << "Tokenized input file list was empty!" << endl;
		print_usage();
		return 1;
	}

	if(optind >= argc) {
		print_usage();
		return 1;
	}
	outfile = argv[optind++];
	vector<String<Dna, Packed<> > > ss;
	
	if(pack) {
			
		try {
			for(size_t i = 0; i < infiles.size(); i++) {
				ifstream in(infiles[i].c_str());
				//TODO: check open success?
				cerr << "packing " << infiles[i]<<"..."<<endl;
				readAndPackFasta(in, ss, verbose);
				in.close();
			}
		} catch(MalformedFastaException& e) {
			cerr << "MalformedFastaException for \"" << infile << "\": " << e.what() << endl;
			return 2;
		} catch(exception& e)
		{
			cerr << "Error packing fasta file \""<< infile << "\": " << e.what() << endl;
			return 2;
		}
		if(verbose) cout << "Read " << ss.size() << " packed sequences" << endl;
		try {
			ofstream out(outfile.c_str());
			writePacked(out, ss, verbose);
			out.close();
		} catch(UnexpectedTypeSizeException& e) {
			cerr << "UnexpectedTypeSizeException for \"" << outfile << "\": " << e.what() << endl;
			return 2;
		} 
		catch(ofstream::failure e) {
			cerr << "Could not write to \"" << outfile << "\": " << e.what() << endl;
			return 2;
		}
		catch(exception& e)
		{
			cerr << "Error packing fasta file \""<< infile << "\": " << e.what() << endl;
			return 2;
		}
	} else {
		unpack(infile, ss, &outfile);
	}

	return 0;
}

#endif
