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
  $Id: pipe_edit_environment.h,v 1.1 2008/08/25 16:20:03 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_PIPE_EDIT_ENVIRONMENT_H
#define SEQAN_HEADER_PIPE_EDIT_ENVIRONMENT_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template < typename TDistanceSpec, unsigned STEP_SIZE = 1 >
    struct EditEnvironment;

/**
.Spec.EditEnvironment:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs tuples of the $tupleLen$ consecutive elements of the input stream.
..signature:Pipe<TInput, Tupler<tupleLen, omitLast> >
..param.TInput:The type of the pipeline module this module reads from.
..param.tupleLen:The tuple length.
...remarks:The tuples contain elements $in[i]in[i+1]...in[i+(tupleLen-1)]$.
..param.omitLast:Omit half filled tuples.
..param.omitLast:If $true$, the output stream is $tupleLen-1$ elements shorter than the input stream.
..param.omitLast:If $false$, the lengths are identical and the last tuples are filled with blanks (default constructed elements) for undefined entries.
..remarks:The output type is a @Class.Tuple@ of input elements and length $tupleLen$ (i.e. $Tuple<Value<TInput>::Type, tupleLen>$).
..remarks:The tuples are sequences of the form $in[i]in[i-1]in[i-2]..in[i-tupleLen+1]$. For $omitLast=false$ $i$ begins with 0 and for $omitLast=true$ $i$ begins with $tupleLen-1$.
*/

    //////////////////////////////////////////////////////////////////////////////
    // pipe to enumerate the hamming 1-environment
    template < typename TInput, unsigned STEP_SIZE >
    struct Pipe< TInput, EditEnvironment< Tag<_HammingDistance>, STEP_SIZE > >
    {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

        TInput                      &in;
        typename Value<Pipe>::Type	tmp, orig;
		unsigned					errorPos;		// position of substitution
		unsigned					character;		// replacement character 
		unsigned					skipChar;		// skip the original character 
        
        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			do {
				if (++character < ValueSize<TValue>::VALUE)
					// next replacement value
					assignValueAt(tmp.i2, errorPos, (TValue) character);
				else {
					// next error position					
					assignValueAt(tmp.i2, errorPos, orig.i2[errorPos]);
					character = 0;
					if (++errorPos < length(tmp.i2)) {
						skipChar = (unsigned) orig.i2[errorPos];
						assignValueAt(tmp.i2, errorPos, (TValue) 0);
					} else {
						// next tuple
						errorPos = 0;
						++in;
						for(unsigned i = 1; i < STEP_SIZE && !eof(in); ++i)
							++in;
						if (!eof(in)) {
							tmp = orig = *in;
							assignValueAt(tmp.i2, 0, (TValue) 0);
						}
					}
				}
			// output the original tuple only once
			} while ((errorPos > 0) && (character == skipChar)); 

            return *this;
        }
	};


    //////////////////////////////////////////////////////////////////////////////
    // pipe to enumerate the levenshtein 1-environment
    template < typename TInput, unsigned STEP_SIZE >
    struct Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance>, STEP_SIZE > >
    {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type	TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		enum TState { _SUBST, _DELETE, _INSERT, _INSERT_LAST, _EOF, _INSERT_EOS };

        TInput                      &in;
        typename Value<Pipe>::Type	tmp, orig, prev;
		unsigned					errorPos;		// position of substitution
		unsigned					character;		// replacement character 
		unsigned					skipChar;		// skip the original character 
		TState						state;

        Pipe(TInput& _in):
            in(_in),
			state(_EOF) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
			switch (state) {
			case _SUBST:
				// before _SUBST (tmp[1..] == orig[1..] and tmp[0] == 0) holds
				do {
					if (++character < ValueSize<TValue>::VALUE) {
						// next replacement value
						assignValueAt(tmp.i2, errorPos, (TValue) character);
					} else {
						// next substitution position
						assignValueAt(tmp.i2, errorPos, orig.i2[errorPos]);
						character = 0;
						if (++errorPos < length(tmp.i2)) {
							skipChar = (unsigned) orig.i2[errorPos];
							assignValueAt(tmp.i2, errorPos, (TValue) 0);
						} else {
							// NEXT TUPLE
							// now (tmp == orig) holds
							++in;
							for(unsigned i = 1; i < STEP_SIZE && !eof(in); ++i)
								++in;
							if (!eos(in)) {
								prev = orig;
								orig = *in;
								tmp.i2 = orig.i2;
								assignValueAt(tmp.i2, 0, prev.i2[0]);
								assignValueAt(tmp.i2, 1, prev.i2[1]);
								if (length(tmp.i2) >= 4) {
									errorPos = 2;
									state = _DELETE;
									//::std::cerr << ::std::endl << "_DELETIONS____" << ::std::endl;
									return *this;
								}
							} else {
								// LAST TUPLE
								shiftLeft(orig.i2);
								assignValueAt(tmp.i2, 0, orig.i2[0]);
								assignValueAt(tmp.i2, 1, (TValue) 0);
								character = 0;
								errorPos = 1;
								state = _INSERT_LAST;
								//::std::cerr << ::std::endl << "_INSERTS______" << ::std::endl;
								return *this;
							}
						}
					}
				// output the original tuple only once
				} while ((errorPos > 0) && (character == skipChar));
				break;
			case _DELETE:
				// before _DELETE (prev=orig, ++in; tmp=orig=*in) holds
				assignValueAt(tmp.i2, errorPos, prev.i2[errorPos]);
				if (++errorPos >= length(tmp.i2) - 1) {
					assignValueAt(tmp.i2, length(tmp.i2)-1, prev.i2[length(tmp.i2)-1]);
					assignValueAt(tmp.i2, 0, orig.i2[0]);
					assignValueAt(tmp.i2, 1, (TValue) 0);
					character = 0;
					errorPos = 1;
					state = _INSERT;
					//::std::cerr << ::std::endl << "_INSERTS______" << ::std::endl;
				}
				break;

			case _INSERT_EOS:
				state = _INSERT;
			case _INSERT:
			case _INSERT_LAST:
				// before _INSERT (prev=orig, ++in; tmp=prev) holds
				if (++character < ValueSize<TValue>::VALUE)
					// next replacement value
					assignValueAt(tmp.i2, errorPos, (TValue) character);
				else {
					// next insert position					
					assignValueAt(tmp.i2, errorPos, orig.i2[errorPos]);
					character = 0;
					if (++errorPos >= length(tmp.i2) - 1 && state == _INSERT) {
						tmp = orig;
						state = _SUBST;
						//::std::cerr << ::std::endl << "_REPLACEMENTS_" << ::std::endl;
						errorPos = 0;
						assignValueAt(tmp.i2, 0, (TValue) 0);
						break;
					}
					if (errorPos >= length(tmp.i2)) {
						if (eof(in))
							state = _EOF;
						else {
							tmp = orig = *in;

							// begin to insert the first char at position 0
							shiftRight(tmp.i2);
							assignValueAt(tmp.i2, 1, (TValue) 0);
							errorPos = 0;
							state = _INSERT_EOS;
							//::std::cerr << ::std::endl << "_INSERTS______" << ::std::endl;
						}
						break;
					}
					assignValueAt(tmp.i2, errorPos, (TValue) 0);
				}
			default:;
			}			
            return *this;
		}
	};


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned STEP_SIZE >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_HammingDistance>, STEP_SIZE > > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;

		me.tmp = me.orig = *me.in;
		me.errorPos = 0;
		me.character = 0;

		return true;
	}
    
    template < typename TInput, unsigned STEP_SIZE >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance>, STEP_SIZE > > &me, 
		ControlBeginRead const &command) 
	{
        if (!control(me.in, command)) return false;

		if (eof(me.in)) {
			me.state = me._EOF;
			return true;
		}

		typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
		typedef typename Value<TTuple>::Type							TValue;

		me.tmp = me.orig = *me.in;

		// begin to insert the first char at position 0
		shiftRight(me.tmp.i2);
		assignValueAt(me.tmp.i2, 0, (TValue) 0);
		me.character = 0;
		me.errorPos = 0;
		me.state = me._INSERT;
		//::std::cerr << ::std::endl << "_INSERTS______" << ::std::endl;

		return true;
	}
    
    template < typename TInput, unsigned STEP_SIZE >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance>, STEP_SIZE > > &me, 
		ControlEof const &) 
	{
		return me.state == me._EOF;
    }

    template < typename TInput, unsigned STEP_SIZE >
	inline bool 
	control(
		Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance>, STEP_SIZE > > &me, 
		ControlEos const &) 
	{
		return me.state == me._EOF || me.state == me.INSERT_EOS;
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<_HammingDistance>, STEP_SIZE > > > >::Type
    length(Pipe< TInput, EditEnvironment< Tag<_HammingDistance>, STEP_SIZE > > const &me) {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
		typedef typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<_HammingDistance>, STEP_SIZE > > > >::Type TSize;

		unsigned alphabetSize = ValueSize< typename Value<TTuple>::Type >::VALUE;
		unsigned seqs = countSequences(me.in);
		TSize sum = 0;
/*		for(unsigned i = 0; i < seqs; ++i)
			sum += (length((*me.in.in.in.in.set)[i]) / STEP_SIZE) * (1 + length(me.tmp.i2) * (alphabetSize - 1));
*/		return (length(me.in) / STEP_SIZE) * (1 + length(me.tmp.i2) * (alphabetSize - 1));
		return sum;
    }

    template < typename TInput, unsigned STEP_SIZE >
    inline typename Size< Pipe< TInput, Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance>, STEP_SIZE > > > >::Type
    length(Pipe< TInput, EditEnvironment< Tag<_LevenshteinDistance>, STEP_SIZE > > const &me) {
		typedef typename Value< typename Value<TInput>::Type, 2 >::Type TTuple;
		unsigned alphabetSize = ValueSize< typename Value<TTuple>::Type >::VALUE;
		unsigned seqs = countSequences(me.in);
		// TODO: We run into problems when one sequence contains 1 or less tuples
		// length should be ommitted in future, but Pools or the skew algorithm needs to know the stream length
		if (length(me.in) >= seqs)
			return 
				  (length(me.in) / STEP_SIZE)     * (1 + length(me.tmp.i2) * (alphabetSize - 1)) +			// substitutions and original
				 ((length(me.in) / STEP_SIZE) - seqs) * (length(me.tmp.i2) - 3) +							// deletions
				(((length(me.in) / STEP_SIZE) + seqs) * (length(me.tmp.i2) - 2) + 2 * seqs) * alphabetSize;	// insertions
		else
			return 0;
    }
//}

}

#endif
