/*
 * ebwt.cpp
 *
 *  Created on: Sep 23, 2009
 *      Author: Ben Langmead
 */

#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

/**
 * Try to find the Bowtie index specified by the user.  First try the
 * exact path given by the user.  Then try the user-provided string
 * appended onto the path of the "indexes" subdirectory below this
 * executable, then try the provided string appended onto
 * "$BOWTIE_INDEXES/".
 */
string adjustEbwtBase(const string& cmdline,
					  const string& ebwtFileBase,
					  bool verbose = false)
{
	string str = ebwtFileBase;
	ifstream in;
	if(verbose) cout << "Trying " << str << endl;
	in.open((str + ".1.ebwt").c_str(), ios_base::in | ios::binary);
	if(!in.is_open()) {
		if(verbose) cout << "  didn't work" << endl;
		in.close();
		str = cmdline;
		size_t st = str.find_last_of("/\\");
		if(st != string::npos) {
			str.erase(st);
			str += "/indexes/";
		} else {
			str = "indexes/";
		}
		str += ebwtFileBase;
		if(verbose) cout << "Trying " << str << endl;
		in.open((str + ".1.ebwt").c_str(), ios_base::in | ios::binary);
		if(!in.is_open()) {
			if(verbose) cout << "  didn't work" << endl;
			in.close();
			if(getenv("BOWTIE_INDEXES") != NULL) {
				str = string(getenv("BOWTIE_INDEXES")) + "/" + ebwtFileBase;
				if(verbose) cout << "Trying " << str << endl;
				in.open((str + ".1.ebwt").c_str(), ios_base::in | ios::binary);
				if(!in.is_open()) {
					if(verbose) cout << "  didn't work" << endl;
					in.close();
				} else {
					if(verbose) cout << "  worked" << endl;
				}
			}
		}
	}
	if(!in.is_open()) {
		cerr << "Could not locate a Bowtie index corresponding to basename \"" << ebwtFileBase << "\"" << endl;
		throw 1;
	}
	return str;
}
