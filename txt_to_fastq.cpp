#ifdef TXT_TO_FASTQ_MAIN

/*
 *  txt_to_fastq.cpp
 *  CSAMapper
 *
 *  Created by Cole Trapnell on 7/12/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "sequence_io.h"
#include <seqan/sequence.h>
#include "pat.h"

using namespace std;
using namespace seqan;

int main(int argc, char** argv)
{
	vector<string> queries;
	for (int i = 1; i < argc; ++i)
	{
		queries.push_back(*(argv + i));
	}
//	
	String<Dna, Alloc<> > tmp;
	
	FastqPatternSource<String<Dna, Alloc<> > > patsrc(queries, false, false, NULL, 0, 0, true);
	
	while(patsrc.hasMorePatterns())
	{
		String<Dna, Alloc<> > pat;
		string qual;
		string name;
		patsrc.nextPattern(pat, qual, name);
		cout << "@"  << name << endl;
		cout << pat  << endl;
		cout << "+"  << endl;
		cout << qual << endl;
	}
	return 0;
}

#endif