/*
 * edit.h
 *
 *  Created on: Jul 31, 2009
 *      Author: Ben Langmead
 */

#ifndef EDIT_H_
#define EDIT_H_

#include <iostream>
#include <stdint.h>
#include "assert_helpers.h"
#include "filebuf.h"

/**
 * 3 types of edits; mismatch (substitution), insertion in the
 * reference, deletion in the reference.
 */
enum {
	EDIT_TYPE_MM = 1,
	EDIT_TYPE_SNP,
	EDIT_TYPE_INS,
	EDIT_TYPE_DEL
};

/**
 * Encapsulates an edit between the read sequence and the reference
 * sequence.
 */
struct Edit {

	Edit() : pos(1023) { }

	Edit(int po, int ch, int ty = EDIT_TYPE_MM) :
		chr(ch), qchr(0), type(ty), pos(po) { }

	/**
	 * Write Edit to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
		assert_eq(4, sizeof(Edit));
		fb.writeChars((const char*)this, 4);
	}

	/**
	 * Read Edit from a FileBuf.
	 */
	void deserialize(FileBuf& fb) {
		fb.get((char*)this, 4);
	}

	/**
	 * Edit less-than overload.
	 */
	int operator< (const Edit &rhs) const {
		if(pos < rhs.pos) return 1;
		if(pos > rhs.pos) return 0;
		if(chr < rhs.chr) return 1;
		if(chr > rhs.chr) return 0;
		return (qchr < rhs.qchr)? 1 : 0;
	}

	/**
	 * Edit equals overload.
	 */
	int operator== (const Edit &rhs) const {
		return(pos  == rhs.pos &&
			   chr  == rhs.chr &&
			   qchr == rhs.qchr &&
			   type == rhs.type);
	}

	/**
	 * Return true iff this Edit is initialized.
	 */
	bool initialized() const {
		return pos != 1023;
	}

	uint32_t chr       :  8; // reference character involved (for subst and ins)
	uint32_t qchr      :  8; // read character involved (for subst and del
	uint32_t type      :  4; // 1 -> mm, 2 -> SNP, 3 -> ins, 4 -> del
	uint32_t pos       : 10; // position w/r/t search root
	uint32_t reserved  :  2; // padding

	friend std::ostream& operator<< (std::ostream& os, const Edit& e);
};

#endif /* EDIT_H_ */
