/*
 * edit.h
 *
 *  Created on: Jul 31, 2009
 *      Author: langmead
 */

#ifndef EDIT_H_
#define EDIT_H_

#include <stdint.h>
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
	 * Return true iff this Edit is initialized.
	 */
	bool initialized() const {
		return type != 0;
	}

	uint16_t type      :  2; // 1 -> subst, 2 -> ins, 3 -> del, 0 -> empty
	uint16_t pos       : 10; // position w/r/t search root
	uint16_t chr       :  2; // character involved (for subst and ins)
	uint16_t reserved  :  2; // reserved
};

#endif /* EDIT_H_ */
