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
										unsigned int allowed_differences = default_allowed_diffs) : 
	_text_strs(ss),
	_mer(left_mer_length), 
	_allow_indels(allow_indels),
	_allowed_diffs(allowed_differences),
	_mer_search_params(_mer_hit_sink,
					   _stats,
					   MHP_CHASE_SAMPLE,
					   os)
	{
	}
	
	virtual void search(Ebwt<TStr>& ebwt,
						EbwtSearchStats<TStr>& stats,
						EbwtSearchParams<TStr>& params,
						const TStr& pat, 
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
		
		EbwtSearchState<TStr> s(ebwt, mer, _mer_search_params, 0);
		ebwt.search(s, _mer_search_params);
		
		int max_diffs = _allowed_diffs;
		
		vector<Hit>& mer_hits = _mer_hit_sink.retainedHits();
		for (vector<Hit>::iterator itr = mer_hits.begin(); 
			 itr != mer_hits.end(); 
			 ++itr)
		{
			assert (itr->h.first < _text_strs.size());
			const String<Dna, Packed<> >& ref = _text_strs[itr->h.first];
			
			int ref_hit_start = itr->h.second;
			int ref_hit_end = ref_hit_start + _mer;
			
			if (!revcomp)
			{
				int ref_right_end = ref_hit_end + length(pat) - _mer;
				if (_allow_indels)
					ref_right_end += max_diffs;
				ref_right_end = min((int)length(ref) - 1, ref_right_end);
				
				String<Dna, Packed<> > ref_right_str;
				DnaWord ref_right(0ull, 0);
				if (ref_hit_end < ref_right_end)
				{
					assert(length(ref)/16 <= length(host(ref)));
					ref_right_str = infix(ref, ref_hit_end, ref_right_end);
					ref_right = pack_dna_string(ref_right_str, false);
				}
				
				String<Dna, Packed<> > pat_str = suffix(pat, _mer);
				DnaWord pat_right = pack_dna_string(pat_str, false);
				
				//cerr << "comparing: " << endl;
				//cerr << "\t" << ref_right_str << endl;
				//cerr << "\t" << pat_str << endl;
				
				int diffs;
				int ref_chars_remaining;
				
				if (_allow_indels)
				{
					edit_distance(ref_right, 
								  pat_right, 
								  max_diffs, 
								  false,
								  &diffs,
								  &ref_chars_remaining);
					
					if (diffs <= max_diffs)
					{
						Hit hit = *itr;
						// FIXME: need a clear range here
						hit_sink.reportHit(hit.h, 
										   hit.patId, 
										   pat, 
										   hit.fw, 
										   hit.mms, 
										   hit.oms);
						if (params.multiHitPolicy() == MHP_PICK_1_RANDOM)
							break;
					}
				}
				else
				{
					bitset<max_read_bp> mism = mismatching_bases(ref_right, 
																pat_right, 
																false);
					if (mism.count() <= _allowed_diffs)
					{
						Hit hit = *itr;
						hit_sink.reportHit(hit.h, 
										   hit.patId, 
										   pat, 
										   hit.fw, 
										   hit.mms | (mism << _mer),
										   hit.oms);
						if (params.multiHitPolicy() == MHP_PICK_1_RANDOM)
							break;	
					}
				}
			}
			else
			{
				int mer_start = length(pat) - _mer;
				int ref_left_start = ref_hit_start - mer_start;
				if (_allow_indels)
					mer_start -= max_diffs;
				
				ref_left_start = max(ref_left_start, 0);
				
				
				String<Dna, Packed<> > ref_left_str;
				DnaWord ref_left(0ull, 0);
				if (ref_left_start != ref_hit_start)
				{
					ref_left_str = infix(ref, ref_left_start, ref_hit_start);
					ref_left = pack_dna_string(ref_left_str, true);
				}
				
				String<Dna, Packed<> > pat_str = prefix(pat, mer_start);
				DnaWord pat_left = pack_dna_string(pat_str, true);
				
//								cerr << "comparing: " << endl;
//								cerr << "\t" << ref_left_str << endl;
//								cerr << "\t" << pat_str << endl;
				
				if (_allow_indels)
				{
				
					int diffs;
					int ref_chars_remaining;
					edit_distance(ref_left, 
								  pat_left, 
								  max_diffs, true,
								  &diffs,
								  &ref_chars_remaining);
					if (diffs <= max_diffs)
					{
						Hit hit = *itr;
						pair<uint32_t, uint32_t> hit_coords;
						hit_coords = make_pair<uint32_t, uint32_t>(hit.h.first, ref_left_start + ref_chars_remaining);
						
						hit_sink.reportHit(hit_coords,
										   hit.patId, 
										   pat,
										   hit.fw,
										   hit.mms, 
										   hit.oms);
						if (params.multiHitPolicy() == MHP_PICK_1_RANDOM)
							break;
					}
				}
				else
				{
					
					bitset<max_read_bp> mism = mismatching_bases(ref_left, 
																pat_left, 
																true);
					if (mism.count() <= _allowed_diffs)
					{
						Hit hit = *itr;
						pair<uint32_t, uint32_t> hit_coords;
						hit_coords = make_pair<uint32_t, uint32_t>(hit.h.first, ref_left_start);
						
						hit_sink.reportHit(hit_coords,
										   hit.patId, 
										   pat,
										   hit.fw,
										   hit.mms | (mism << _mer), 
										   hit.oms);
						if (params.multiHitPolicy() == MHP_PICK_1_RANDOM)
							break;	
					}
				}
			}
			
		}
		
		_mer_hit_sink.clearRetainedHits();
	};
protected:
	vector<String<Dna, Packed<> > >& _text_strs;
	unsigned int _mer;
	bool _allow_indels;
	unsigned int _allowed_diffs;
	HitBucket _mer_hit_sink;
	EbwtSearchStats<TStr> _stats;
	EbwtSearchParams<TStr> _mer_search_params;
	vector<TStr> os;
};

#endif

