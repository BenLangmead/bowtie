#ifndef TOKENIZE_H_
#define TOKENIZE_H_

#include <string>
#include <vector>

using namespace std;

/**
 * Split string s according to given delimiters.  Mostly borrowed
 * from C++ Programming HOWTO 7.3.
 */
static inline void tokenize(const string& s, const string& delims, vector<string>& ss) {
	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos = s.find_first_of(delims, lastPos);
	while (string::npos != pos || string::npos != lastPos) {
		ss.push_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delims, pos);
        pos = s.find_first_of(delims, lastPos);
    }
}

#endif /*TOKENIZE_H_*/
