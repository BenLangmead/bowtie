/*
 *  hit.cpp
 *  CSAMapper
 *
 *  Created by Cole Trapnell on 6/18/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include "hit.h"

/// Sort by text-id then by text-offset
bool operator< (const Hit& a, const Hit& b) {
    return a.h < b.h;
}
