/*
 * hit_set.h
 *
 *  Created on: Jul 31, 2009
 *      Author: Ben Langmead
 */

#ifndef HIT_SET_H_
#define HIT_SET_H_

#include <seqan/sequence.h>
#include <vector>
#include <algorithm>
#include "assert_helpers.h"
#include "filebuf.h"
#include "edit.h"
#include "alphabet.h"
#include "annot.h"
#include "refmap.h"

/**
 * Encapsulates a hit contained within a HitSet that can be
 * (de)serialized to/from FileBufs.  Used for chaining.
 */
struct HitSetEnt {
	typedef std::pair<uint32_t,uint32_t> U32Pair;

	HitSetEnt() { }

	/**
	 * Write binary representation of HitSetEnt to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
		fb.writeChars((const char*)&h.first, 4);
		fb.writeChars((const char*)&h.second, 4);
		assert(fw == 0 || fw == 1);
		fb.write(fw);
		assert_geq(stratum, 0);
		assert_lt(stratum, 4);
		fb.write(stratum);
		assert_eq(stratum, (cost >> 14));
		fb.writeChars((const char*)&cost, 2);
		fb.writeChars((const char*)&oms, 4);
		uint32_t sz = edits.size();
		fb.writeChars((const char*)&sz, 4);
		std::vector<Edit>::const_iterator it;
		for(it = edits.begin(); it != edits.end(); it++) {
			it->serialize(fb);
		}
		sz = cedits.size();
		fb.writeChars((const char*)&sz, 4);
		for(it = cedits.begin(); it != cedits.end(); it++) {
			it->serialize(fb);
		}
	}

	/**
	 * Repopulate a HitSetEnt from its binary representation in FileBuf.
	 */
	void deserialize(FileBuf& fb) {
		fb.get((char*)&h.first, 4);
		fb.get((char*)&h.second, 4);
		fw = fb.get();
		assert(fw == 0 || fw == 1);
		stratum = fb.get();
		assert_geq(stratum, 0);
		assert_lt(stratum, 4);
		fb.get((char*)&cost, 2);
		assert_eq(stratum, (cost >> 14));
		fb.get((char*)&oms, 4);
		uint32_t sz = 0;
		fb.get((char*)&sz, 4);
		assert_lt(sz, 1024);
		edits.resize(sz);
		for(uint32_t i = 0; i < sz; i++) {
			edits[i].deserialize(fb);
		}
		fb.get((char*)&sz, 4);
		assert_lt(sz, 1024);
		cedits.resize(sz);
		for(uint32_t i = 0; i < sz; i++) {
			cedits[i].deserialize(fb);
		}
	}

	/**
	 * Less than operator.  Break HitSetEnt ties by:
	 *  - Stratum, then
	 *  - Cost, then
	 *  - Position, then
	 *  - Orientation
	 */
	int operator< (const HitSetEnt &rhs) const {
		if(stratum < rhs.stratum) return 1;
		if(stratum > rhs.stratum) return 0;
		if(cost < rhs.cost) return 1;
		if(cost > rhs.cost) return 0;
		if(h < rhs.h) return 1;
		if(h > rhs.h) return 0;
		return (fw < rhs.fw)? 1 : 0;
	}

	/**
	 * Greater than operator.
	 */
	int operator> (const HitSetEnt &rhs) const {
		if(stratum < rhs.stratum) return 0;
		if(stratum > rhs.stratum) return 1;
		if(cost < rhs.cost) return 0;
		if(cost > rhs.cost) return 1;
		if(h < rhs.h) return 0;
		if(h > rhs.h) return 1;
		return (fw <= rhs.fw)? 0 : 1;
	}

	/**
	 * Equality comparison operator.
	 */
	int operator== (const HitSetEnt &rhs) const {
		return(stratum == rhs.stratum &&
		       cost == rhs.cost &&
		       fw == rhs.fw &&
		       h == rhs.h);
	}

	/**
	 * Indexing returns edits.
	 */
	Edit& operator[](unsigned x) {
		return edits[x];
	}

	/**
	 * Indexing returns edits.
	 */
	const Edit& operator[](unsigned x) const {
		return edits[x];
	}

	/**
	 * Another way to get at an edit.
	 */
	Edit& editAt(unsigned i) {
		return edits[i];
	}

	/**
	 * Another way to get at a const edit.
	 */
	const Edit& editAt(unsigned i) const {
		return edits[i];
	}

	/**
	 * Get the ith color edit.
	 */
	Edit& colorEditAt(unsigned i) {
		return cedits[i];
	}

	/**
	 * Another way to get at an edit.
	 */
	const Edit& colorEditAt(unsigned i) const {
		return cedits[i];
	}

	/**
	 * Return the front entry.
	 */
	Edit& front() {
		return edits.front();
	}

	/**
	 * Return the back entry.
	 */
	Edit& back() {
		return edits.back();
	}

	/**
	 * Expand the entry list by one.
	 */
	void expand() {
		edits.resize(edits.size() + 1);
	}

	/**
	 * Sort edits by position
	 */
	void sort() {
		if(edits.size() > 1) std::sort(edits.begin(), edits.end());
	}

	/**
	 * Return number of edits.
	 */
	size_t size() const {
		return edits.size();
	}

	bool empty() const {
		return edits.empty();
	}

	/**
	 * Write HitSetEnt to an output stream.
	 */
	friend std::ostream& operator << (std::ostream& os, const HitSetEnt& hse);

	U32Pair h; // reference coordinates
	uint8_t fw; // orientation
	int8_t stratum; // stratum
	uint16_t cost; // cost, including stratum
	uint32_t oms; // # others
	std::vector<Edit> edits; // edits to get from reference to subject
	std::vector<Edit> cedits; // color edits to get from reference to subject
};

/**
 * Encapsulates a set of hits that can be (de)serialized to/from
 * FileBufs.  Used for chaining.
 */
struct HitSet {

	typedef std::vector<HitSetEnt> EntVec;
	typedef EntVec::const_iterator Iter;
	typedef std::pair<uint32_t,uint32_t> U32Pair;

	HitSet() {
		maxedStratum = -1;
	}

	HitSet(FileBuf& fb) {
		deserialize(fb);
	}

	/**
	 * Write binary representation of HitSet to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
		fb.write(color ? 1 : 0);
		uint32_t i = seqan::length(name);
		assert_gt(i, 0);
		fb.writeChars((const char*)&i, 4);
		fb.writeChars(seqan::begin(name), i);
		i = seqan::length(seq);
		assert_gt(i, 0);
		assert_lt(i, 1024);
		fb.writeChars((const char*)&i, 4);
		for(size_t j = 0; j < i; j++) {
			fb.write("ACGTN"[(int)seq[j]]);
		}
		fb.writeChars(seqan::begin(qual), i);
		i = ents.size();
		fb.writeChars((const char*)&i, 4);
		std::vector<HitSetEnt>::const_iterator it;
		for(it = ents.begin(); it != ents.end(); it++) {
			it->serialize(fb);
		}
		fb.write(maxedStratum);
	}

	/**
	 * Repopulate a HitSet from its binary representation in FileBuf.
	 */
	void deserialize(FileBuf& fb) {
		color = (fb.get() != 0 ? true : false);
		uint32_t sz = 0;
		if(fb.get((char*)&sz, 4) != 4) {
			seqan::clear(name);
			seqan::clear(seq);
			return;
		}
		assert_gt(sz, 0);
		assert_lt(sz, 1024);
		seqan::resize(name, sz);
		fb.get(seqan::begin(name), sz);
		fb.get((char*)&sz, 4);
		assert_gt(sz, 0);
		assert_lt(sz, 1024);
		seqan::resize(seq, sz);
		for(size_t j = 0; j < sz; j++) {
			seq[j] = charToDna5[fb.get()];
		}
		seqan::resize(qual, sz);
		fb.get(seqan::begin(qual), sz);
		fb.get((char*)&sz, 4);
		if(sz > 0) {
			ents.resize(sz);
			for(size_t i = 0; i < sz; i++) {
				ents[i].deserialize(fb);
			}
		} else {
			ents.clear();
		}
		maxedStratum = fb.get();
	}

	/**
	 * Return true iff this HitSet is initialized with a read (but not
	 * necessarily any alignments).
	 */
	bool initialized() const {
		return !seqan::empty(seq);
	}

	/**
	 * Return true iff this HitSet has no hits.
	 */
	bool empty() const {
		return ents.empty();
	}

	/**
	 * Return number of entries in this HitSet.
	 */
	size_t size() const {
		return ents.size();
	}

	/**
	 * Remove all entries from this HitSet.
	 */
	void clear() {
		maxedStratum = -1;
		ents.clear();
	}

	Iter begin() const {
		return ents.begin();
	}

	Iter end() const {
		return ents.end();
	}

	/**
	 * Return the front entry.
	 */
	HitSetEnt& front() {
		return ents.front();
	}

	/**
	 * Return the back entry.
	 */
	HitSetEnt& back() {
		return ents.back();
	}

	/**
	 * Expand the entry list by one.
	 */
	void expand() {
		ents.resize(ents.size() + 1);
		assert(ents.back().empty());
	}

	/**
	 * Resize entry list
	 */
	void resize(size_t sz) {
		ents.resize(sz);
	}

	/**
	 * Sort hits by stratum/penalty.
	 */
	void sort() {
		if(ents.size() > 1) std::sort(ents.begin(), ents.end());
	}

	/**
	 * Return true iff Ents are sorted
	 */
	bool sorted() const {
		if(ents.empty()) return true;
		for(size_t i = 0; i < ents.size()-1; i++) {
			if(!(ents[i] < ents[i+1])) return false;
		}
		return true;
	}

	/**
	 * Remove the ith hit from the HitSet, shifting all subsequent hits
	 * up by one.
	 */
	void remove(size_t i) {
		ents.erase(ents.begin() + i);
	}

	/**
	 * Return true if strata are uniform across hits; assert if they're
	 * not.
	 */
	bool uniformStrata() const {
		for(size_t i = 0; i < ents.size()-1; i++) {
			assert(ents[i].stratum == ents[i+1].stratum);
		}
		return true;
	}

	/**
	 * Indexing returns entries.
	 */
	HitSetEnt& operator[](unsigned x) {
		return ents[x];
	}

	/**
	 * Indexing returns entries.
	 */
	const HitSetEnt& operator[](unsigned x) const {
		return ents[x];
	}

	/**
	 * Apply a reference mappings to all the contained hits.
	 */
	void applyReferenceMap(const ReferenceMap& map) {
		std::vector<HitSetEnt>::iterator it;
		for(it = ents.begin(); it != ents.end(); it++) map.map(it->h);
	}

	/**
	 * Clear out all the strings and all the entries.
	 */
	void clearAll() {
		seqan::clear(name);
		seqan::clear(seq);
		seqan::clear(qual);
		ents.clear();
		color = false;
	}

	/**
	 *
	 */
	bool tryReplacing(U32Pair h,
	                  bool fw,
	                  uint16_t cost,
	                  size_t& pos)
	{
		for(size_t i = 0; i < ents.size(); i++) {
			if(ents[i].h == h && ents[i].fw == fw) {
				if(cost < ents[i].cost) {
					// New hit at same position is better!  Replace it.
					ents[i].h = h;
					ents[i].fw = fw;
					ents[i].stratum = cost >> 14;
					ents[i].cost = cost;
					pos = i;
					return true;
				} else {
					pos = 0xffffffff;
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Report up to 'khits' hits from this HitSet.
	 */
	void reportUpTo(std::ostream& os, int khits = INT_MAX);

	/**
	 * Print this HitSet.
	 */
	friend std::ostream& operator << (std::ostream& os, const HitSet& hs);

	seqan::String<char> name;
	seqan::String<seqan::Dna5> seq;
	seqan::String<char> qual;
	int8_t maxedStratum;
	std::vector<HitSetEnt> ents;
	bool color; // whether read was orginally in colorspace
};

#endif /* HIT_SET_H_ */
