/*
 * edit.cpp
 *
 *  Created on: Jul 14, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include "edit.h"

using namespace std;

ostream& operator<< (ostream& os, const Edit& e) {
	os << e.pos << "ACGT"[e.chr];
	return os;
}
