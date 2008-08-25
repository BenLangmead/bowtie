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
  $Id: system_base.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SYSTEM_BASE_H
#define SEQAN_HEADER_SYSTEM_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

#ifdef SEQAN_DEBUG

#define SEQAN_DO_SYS(_cond) if (!(_cond)) ::SEQAN_NAMESPACE_MAIN::debug::Message< ::SEQAN_NAMESPACE_MAIN::debug::Check >(__FILE__, __LINE__, #_cond " is FALSE");
#define SEQAN_DO_SYS1(_cond) SEQAN_DO_SYS(_cond)
#define SEQAN_DO_SYS2(_cond, _comment) if (!(_cond)) ::SEQAN_NAMESPACE_MAIN::debug::Error< ::SEQAN_NAMESPACE_MAIN::debug::Check >(__FILE__, __LINE__, _comment);

#else

#define SEQAN_DO_SYS(_cond) { (_cond); }
#define SEQAN_DO_SYS1(_cond) SEQAN_DO_SYS(_cond)
#define SEQAN_DO_SYS2(_cond, _comment) SEQAN_DO_SYS(_cond)

#endif

}

#endif
