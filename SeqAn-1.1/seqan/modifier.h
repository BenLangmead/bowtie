 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: modifier.h,v 1.1 2008/08/25 16:20:06 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_H
#define SEQAN_HEADER_MODIFIER_H

//____________________________________________________________________________
// prerequisites

#include <functional>

//____________________________________________________________________________
// basics

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/modifier/modifier_generated_forwards.h>
#endif

#include <seqan/sequence.h>
#include <seqan/modifier/modifier_iterator.h>
#include <seqan/modifier/modifier_string.h>

//____________________________________________________________________________
// applications

#include <seqan/modifier/modifier_functors.h>
#include <seqan/modifier/modifier_view.h>
#include <seqan/modifier/modifier_reverse.h>
#include <seqan/modifier/modifier_shortcuts.h>


#endif //#ifndef SEQAN_HEADER_...
