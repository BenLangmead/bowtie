// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Facade header for module sequence.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_SEQUENCE_H_
#define EXTRAS_INCLUDE_SEQAN_SEQUENCE_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/sequence.h>

// ===========================================================================
// Prerequisites from extras.
// ===========================================================================

#include <seqan/basic_extras.h>

// ===========================================================================
// Container View.
// ===========================================================================

#include <seqan/sequence/container_view.h>

// ===========================================================================
// Adaption of thrust::device_vector.
// ===========================================================================

#ifdef PLATFORM_CUDA
#include <seqan/sequence/adapt_thrust_vector.h>
#endif

// ===========================================================================
// StringSet View and Device.
// ===========================================================================

#include <seqan/sequence/string_set_contat_direct_view.h>
#ifdef PLATFORM_CUDA
#include <seqan/sequence/string_set_contat_direct_device.h>
#endif


#endif  // EXTRAS_INCLUDE_SEQAN_SEQUENCE_H_
