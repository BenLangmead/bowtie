/*
 * tokenize.h
 *
 *  Created on: Jul 21, 2009
 *      Author: Ben Langmead
 */

#ifndef TOKENIZE_H_
#define TOKENIZE_H_

#include <string>
#include <sstream>
#include <vector>
#include <limits>

using namespace std;

/**
 * Split string s according to given delimiters.  Mostly borrowed
 * from C++ Programming HOWTO 7.3.
 */
static inline void tokenize(
	const string& s,
	const string& delims,
	vector<string>& ss,
	size_t max = std::numeric_limits<size_t>::max())
{
	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos = s.find_first_of(delims, lastPos);
	while (string::npos != pos || string::npos != lastPos) {
		ss.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
		if(ss.size() == (max - 1)) {
			pos = string::npos;
		}
	}
}

static inline void tokenize(
	const std::string& s,
    char delim,
    std::vector<std::string>& ss)
{
	std::string token;
	std::istringstream iss(s);
	while(getline(iss, token, delim)) {
		ss.push_back(token);
	}
}

#endif /*TOKENIZE_H_*/
