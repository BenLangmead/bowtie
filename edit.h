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
	EDIT_TYPE_MM = 0,
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
		type(ty), pos(po), chr(ch) { }

	/**
	 * Write Edit to an OutFileBuf.
	 */
	void serialize(OutFileBuf& fb) const {
		assert_eq(2, sizeof(Edit));
		fb.writeChars((const char*)this, 2);
	}

	/**
	 * Read Edit from a FileBuf.
	 */
	void deserialize(FileBuf& fb) {
		fb.get((char*)this, 2);
	}

	/**
	 * Edit less-than overload.
	 */
	int operator< (const Edit &rhs) const {
		if(pos < rhs.pos) return 1;
		if(pos > rhs.pos) return 0;
		return (chr < rhs.chr)? 1 : 0;
	}

	/**
	 * Edit equals overload.
	 */
	int operator== (const Edit &rhs) const {
		return(pos  == rhs.pos &&
			   chr  == rhs.chr &&
			   type == rhs.type);
	}

	/**
	 * Return true iff this Edit is initialized.
	 */
	bool initialized() const {
		return pos != 1023;
	}

	uint16_t type      :  2; // 1 -> subst, 2 -> ins, 3 -> del, 0 -> empty
	uint16_t pos       : 10; // position w/r/t search root
	uint16_t chr       :  2; // character involved (for subst and ins)
	uint16_t reserved  :  2; // reserved

	friend std::ostream& operator<< (std::ostream& os, const Edit& e);
};

#endif /* EDIT_H_ */
