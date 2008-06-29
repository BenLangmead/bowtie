/*
 *  packed_io.cpp
 *  CSAMapper
 *
 *  Created by Cole Trapnell on 6/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "packed_io.h"

void unpack(const string& infile,
			vector<String<Dna, Packed<> > >& ss,
			string* outfile)
{
	// Unpack
	ifstream in(infile.c_str());
	try {
		readPacked(in, ss, -1, false);
	} catch(MalformedFastaException& e) {
		cerr << "MalformedFastaException for \"" << infile << "\": " << e.what() << endl;
		return;
	}
	in.close();
	
	if(outfile)
	{
		ofstream out(outfile->c_str());
		for(size_t i = 0; i < ss.size(); i++) {
			out << ">" << endl;
			out << ss[i] << endl;
		}
		out.close();
	}
} 
