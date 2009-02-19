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
  $Id: pipe.h,v 1.2 2009/02/19 01:51:23 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_H
#define SEQAN_HEADER_PIPE_H

//____________________________________________________________________________
// prerequisites

#include <seqan/file.h>
#include <seqan/basic/basic_volatile_ptr.h>

#include <cstdio>
#include <cassert>
#include <functional>
#include <iterator>
#include <climits>
#include <vector>
#include <queue>

//____________________________________________________________________________
// pipes

#define SEQAN_NAMESPACE_PIPELINING pipe

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/pipe/pipe_generated_forwards.h>
#endif

#include <seqan/pipe/pipe_base.h>
#include <seqan/pipe/pipe_iterator.h>
#include <seqan/pipe/pipe_caster.h>
#include <seqan/pipe/pipe_counter.h>
#include <seqan/pipe/pipe_echoer.h>
#include <seqan/pipe/pipe_edit_environment.h>
#include <seqan/pipe/pipe_filter.h>
#include <seqan/pipe/pipe_joiner.h>
#include <seqan/pipe/pipe_namer.h>
#include <seqan/pipe/pipe_sampler.h>
#include <seqan/pipe/pipe_shifter.h>
#include <seqan/pipe/pipe_source.h>
#include <seqan/pipe/pipe_tupler.h>

//____________________________________________________________________________
// pools

#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

#endif //#ifndef SEQAN_HEADER_...
