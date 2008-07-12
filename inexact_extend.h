#ifndef INEXACT_EXTEND_H
#define INEXACT_EXTEND_H
/*
 *  inexact_extend.h
 *  Bowtie
 *
 *  Created by Cole Trapnell on 6/18/08.
 *  
 *
 */
#include <vector>
#include <set>
#include <utility>
#include <algorithm>
#include <cassert>
#include "hit.h"
#include "ebwt.h"
#include "sequence_io.h"
#include "LVKernel.h"
#include "seqan/sequence.h"

using namespace std;
using std::set;

DnaWord pack_dna_string(String<Dna,Packed<> >& dna, bool pack_left);

template<class TStr>
class SearchPolicy
{
public:
	virtual void search(Ebwt<TStr>& ebwt,
						EbwtSearchStats<TStr>& stats,
						EbwtSearchParams<TStr>& params,
						const TStr& pat,
						const string& patName,
						const string& patQuals,
						HitSink& hit_sink) = 0;
	virtual ~SearchPolicy() {}
};

template<class TStr>
class ExactSearch : public SearchPolicy<TStr>
{
public:
	virtual void search(Ebwt<TStr>& ebwt,
						EbwtSearchStats<TStr>& stats,
						EbwtSearchParams<TStr>& params,
						const TStr& pat,
						const string& patName,
						const string& patQuals,
						HitSink& hit_sink) 
	{
		EbwtSearchState<TStr> s(ebwt, pat, params, 0);
		ebwt.search(s);
	};
};

bitset<max_read_bp> mismatching_bases(const DnaWord& w1, 
									  const DnaWord& w2, 
									  bool left_extend);

static const unsigned int default_allowed_diffs = 4;
static const unsigned int default_left_mer_length = 22;

template<class TStr>
class ExactSearchWithLowQualityThreePrime : public SearchPolicy<TStr>
{
public:
	ExactSearchWithLowQualityThreePrime(vector<String<Dna, Packed<> > >& ss,
										bool allow_indels = false,
										unsigned int left_mer_length = default_left_mer_length,
										unsigned int allowed_differences = default_allowed_diffs,
										bool allow_revcomp = true) : 
	_text_strs(ss),
	_mer(left_mer_length), 
	_allow_indels(allow_indels),
	_allowed_diffs(allowed_differences),
	_mer_search_params(_mer_hit_sink,
					   _stats,
					   MHP_CHASE_SAMPLE,
					   os,
					   allow_revcomp)
	{
	}
	
	virtual void search(Ebwt<TStr>& ebwt,
						EbwtSearchStats<TStr>& stats,
						EbwtSearchParams<TStr>& params,
						const TStr& pat,
						const string& patName,
						const string& patQuals,
						HitSink& hit_sink) 
	{
		
		TStr mer;
		bool revcomp = !params.fw();
		if (!revcomp)
			mer = prefix(pat, _mer);
		else
			mer = suffix(pat, length(pat) - _mer);
		//cerr << mer << endl;
		_mer_search_params.setPatId(params.patId());
		_mer_search_params.setFw(!revcomp);
		_mer_search_params.setEbwtFw(params.ebwtFw());
		
		EbwtSearchState<TStr> s(ebwt, mer, patName, "", _mer_search_params, 0);
		ebwt.search(s, _mer_search_params);
		
		extendHits(pat, patQuals, hit_sink, params, revcomp, revcomp, false);
		
		_mer_hit_sink.clearRetainedHits();
	}
	
	void extendHits(const TStr& pat,
					const string& patQuals,
					HitSink& hit_sink, 
					EbwtSearchParams<TStr>& params,
					bool extend_left_pat,
					bool extend_left_ref,
					bool swap_mismatches)
	{
		TStr pat_to_print = pat;
		
		if (extend_left_ref != extend_left_pat)
			::reverse(pat_to_print);
		
		vector<Hit>& mer_hits = _mer_hit_sink.retainedHits();
		for (vector<Hit>::iterator itr = mer_hits.begin(); 
			 itr != mer_hits.end(); 
			 ++itr)
		{
			assert (itr->h.first < _text_strs.size());
			const String<Dna, Packed<> >& ref = _text_strs[itr->h.first];
			
			int ref_hit_start = itr->h.second;
			int ref_hit_end = ref_hit_start + _mer;
			
			DnaWord pat_word(0ull,0);
			DnaWord ref_word(0ull,0);
			String<Dna, Packed<> > pat_str;
			String<Dna, Packed<> > ref_str;
			if (!extend_left_ref)
			{
				int ref_right_end = ref_hit_end + length(pat) - _mer;
				
				ref_right_end = min((int)length(ref) - 1, ref_right_end);
				
				if (ref_hit_end < ref_right_end)
				{
					assert(length(ref)/16 <= length(host(ref)));
					ref_str = infix(ref, ref_hit_end, ref_right_end);
					ref_word = pack_dna_string(ref_str, false);
				}
			}
			else
			{
				int mer_start = length(pat) - _mer;
				int ref_left_start = ref_hit_start - mer_start;
				if (_allow_indels)
					mer_start -= _allowed_diffs;
				
				ref_left_start = max(ref_left_start, 0);
				
				if (ref_left_start != ref_hit_start)
				{
					ref_str = infix(ref, ref_left_start, ref_hit_start);
					ref_word = pack_dna_string(ref_str, false);
					ref_hit_start = ref_left_start;
				}
			}
				
			if (!extend_left_pat)
			{
				pat_str = suffix(pat, _mer);
				//cerr << "\t len: " <<length(pat)<<" "<< pat_str << endl;
				if (extend_left_ref)
					::reverse(pat_str);
				pat_word = pack_dna_string(pat_str, false);
			}
			else
			{
				pat_str = prefix(pat, length(pat) - _mer);
				//cerr << "\t" << pat_str << endl;
				if (!extend_left_ref)
					::reverse(pat_str);
				pat_word = pack_dna_string(pat_str, false);	
			}
//			cerr << "pattern: "<< pat << endl;
//			cerr << "comparing: " << endl;
//			cerr << "\t" << ref_str << endl;
//			cerr << "\t" << pat_str << endl;
			
			bitset<max_read_bp> mism = mismatching_bases(ref_word, 
														 pat_word, 
														 false);
			
			if (mism.count() <= _allowed_diffs)
			{
				bitset<max_read_bp> mism_to_print;
				unsigned int ext_len;
				
				//if (!swap_mismatches)
//				{
//					ext_len = length(pat);
//					mism_to_print = mism;
//					//cerr << mism_to_print << endl;
//
//					for (unsigned int i = 0; i < mism.size() && i < ext_len; ++i)
//					{
//						if (mism.test(i) && !mism.test(ext_len - i - 1))
//						{
//							mism_to_print.reset(i);
//							mism_to_print.set(ext_len - i - 1,1);
//						}
//					}
//					mism = mism_to_print;
//					//cerr << mism_to_print << endl;
//				}
//				if(extend_left_pat)
//				{
//					mism_to_print = itr->mms;
//					//cerr << mism_to_print << endl;
//					ext_len = _mer;
//					for (unsigned int i = 0; i < itr->mms.size() && i < ext_len; ++i)
//					{
//						if (itr->mms.test(i) && !itr->mms.test(ext_len - i - 1))
//						{
//							mism_to_print.reset(i);
//							mism_to_print.set(ext_len - i - 1,1);
//						}
//					}
//					itr->mms = mism_to_print;
//					//cerr << mism_to_print << endl;
////					if (extend_left_ref)
////					{
////						mism_to_print = itr->mms;
////						//cerr << mism_to_print << endl;
////						ext_len = length(pat);
////						for (unsigned int i = 0; i < itr->mms.size() && i < ext_len; ++i)
////						{
////							if (itr->mms.test(i) && !itr->mms.test(ext_len - i - 1))
////							{
////								mism_to_print.reset(i);
////								mism_to_print.set(ext_len - i - 1,1);
////							}
////						}
////						itr->mms = mism_to_print;
////						//cerr << mism_to_print << endl;
////					}
//				}
					
//					mism_to_print = itr->mms;
//					ext_len = length(pat);
//					for (unsigned int i = 0; i < itr->mms.size() && i < ext_len; ++i)
//					{
//						if (itr->mms.test(i) && !itr->mms.test(ext_len - i - 1))
//						{
//							mism_to_print.reset(i);
//							mism_to_print.set(ext_len - i - 1,1);
//						}
//					}
//					itr->mms = mism_to_print;
					
				
				Hit hit = *itr;
				hit.h.second = ref_hit_start;
				
				bitset<max_read_bp> mms;
				if (!swap_mismatches)
					mms = hit.mms | (mism << _mer);
				else
					mms = mism | (hit.mms << (length(pat) - _mer));
				
				hit_sink.reportHit(hit.h, 
								   hit.patId,
								   hit.patName,
								   pat_to_print, 
								   patQuals,
								   hit.fw, 
								   mms,
								   hit.oms);
				if (params.multiHitPolicy() == MHP_PICK_1_RANDOM)
					break;	
			}
		}
		
	}
	
	vector<String<Dna, Packed<> > >& _text_strs;
	unsigned int _mer;
	bool _allow_indels;
	unsigned int _allowed_diffs;
	HitBucket _mer_hit_sink;
	EbwtSearchStats<TStr> _stats;
	EbwtSearchParams<TStr> _mer_search_params;
	vector<TStr> os;
};

template<class TStr>
class OneMismatchSearchWithLowQualityThreePrime : public ExactSearchWithLowQualityThreePrime<TStr>
{
	
public:
	OneMismatchSearchWithLowQualityThreePrime(vector<String<Dna, Packed<> > >& ss,
											  bool allow_indels = false,
											  unsigned int left_mer_length = default_left_mer_length,
											  unsigned int allowed_differences = default_allowed_diffs,
											  bool allow_revcomp = true) : 
	ExactSearchWithLowQualityThreePrime<TStr>(ss, allow_indels, left_mer_length, allowed_differences, allow_revcomp),
	_suppressExact(false)
	{
	}
	
	void suppressExact(bool _supExact) { _suppressExact = _supExact; } 
	
	virtual void search(Ebwt<TStr>& ebwt,
						EbwtSearchStats<TStr>& stats,
						EbwtSearchParams<TStr>& params,
						const TStr& pat,
						const string& patName,
						const string& patQuals,
						HitSink& hit_sink)
	{
		
		TStr mer;
		//bool revcomp = !params.fw();
		
		bool extend_left_pat = params.ebwtFw() != params.fw();
		//bool extend_left_ref = 
		
		if (extend_left_pat)
			mer = suffix(pat, length(pat) - this->_mer);
		else
			mer = prefix(pat, this->_mer);
		
		//cerr << "searching for: " << mer << " len: " << length(mer) << endl;
		this->_mer_search_params.setPatId(params.patId());
		this->_mer_search_params.setFw(params.fw());
		this->_mer_search_params.setEbwtFw(params.ebwtFw());
		
		
		EbwtSearchState<TStr> s(ebwt, mer, patName, "", this->_mer_search_params, 0);
		ebwt.search1MismatchOrBetter(s, this->_mer_search_params, !_suppressExact);
		
		extendHits(pat, patQuals, hit_sink, params, extend_left_pat, !params.fw(), !params.fw());
		
		this->_mer_hit_sink.clearRetainedHits();
	};
protected:
	bool _suppressExact;
};

#endif

