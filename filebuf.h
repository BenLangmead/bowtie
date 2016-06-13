/*
 * filebuf.h
 *
 *      Author: Ben Langmead
 */
#ifndef FILEBUF_H_
#define FILEBUF_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdexcept>
#include "assert_helpers.h"

/**
 * Simple wrapper for a FILE*, istream or ifstream that reads it in
 * chunks (with fread) and keeps those chunks in a buffer.  It also
 * services calls to get(), peek() and gets() from the buffer, reading
 * in additional chunks when necessary.
 */
class FileBuf {
public:
	FileBuf() {
		init();
	}

	FileBuf(FILE *in) {
		init();
		_in = in;
		assert(_in != NULL);
	}

	FileBuf(std::ifstream *inf) {
		init();
		_inf = inf;
		assert(_inf != NULL);
	}

	FileBuf(std::istream *ins) {
		init();
		_ins = ins;
		assert(_ins != NULL);
	}

	bool isOpen() {
		return _in != NULL || _inf != NULL || _ins != NULL;
	}

	/**
	 * Close the input stream (if that's possible)
	 */
	void close() {
		if(_in != NULL && _in != stdin) {
			fclose(_in);
		} else if(_inf != NULL) {
			_inf->close();
		} else {
			// can't close _ins
		}
	}

	/**
	 * Get the next character of input and advance.
	 */
	int get() {
		assert(_in != NULL || _inf != NULL || _ins != NULL);
		int c = peek();
		if(c != -1) {
			_cur++;
			if(_lastn_cur < LASTN_BUF_SZ) _lastn_buf[_lastn_cur++] = c;
		}
		return c;
	}

	/**
	 * Return true iff all input is exhausted.
	 */
	bool eof() {
		return (_cur == _buf_sz) && _done;
	}

	/**
	 * Initialize the buffer with a new C-style file.
	 */
	void newFile(FILE *in) {
		_in = in;
		_inf = NULL;
		_ins = NULL;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Initialize the buffer with a new ifstream.
	 */
	void newFile(std::ifstream *__inf) {
		_in = NULL;
		_inf = __inf;
		_ins = NULL;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Initialize the buffer with a new istream.
	 */
	void newFile(std::istream *__ins) {
		_in = NULL;
		_inf = NULL;
		_ins = __ins;
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Restore state as though we just started reading the input
	 * stream.
	 */
	void reset() {
		if(_inf != NULL) {
			_inf->clear();
			_inf->seekg(0, std::ios::beg);
		} else if(_ins != NULL) {
			_ins->clear();
			_ins->seekg(0, std::ios::beg);
		} else {
			rewind(_in);
		}
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	/**
	 * Peek at the next character of the input stream without
	 * advancing.  Typically we can simple read it from the buffer.
	 * Occasionally we'll need to read in a new buffer's worth of data.
	 */
	int peek() {
		assert(_in != NULL || _inf != NULL || _ins != NULL);
		assert_leq(_cur, _buf_sz);
		if(_cur == _buf_sz) {
			if(_done) {
				// We already exhausted the input stream
				return -1;
			}
			// Read a new buffer's worth of data
			else {
				// Get the next chunk
				if(_inf != NULL) {
					_inf->read((char*)_buf, BUF_SZ);
					_buf_sz = _inf->gcount();
				} else if(_ins != NULL) {
					_ins->read((char*)_buf, BUF_SZ);
					_buf_sz = _ins->gcount();
				} else {
					assert(_in != NULL);
					_buf_sz = fread(_buf, 1, BUF_SZ, _in);
				}
				_cur = 0;
				if(_buf_sz == 0) {
					// Exhausted, and we have nothing to return to the
					// caller
					_done = true;
					return -1;
				} else if(_buf_sz < BUF_SZ) {
					// Exhausted
					_done = true;
				}
			}
		}
		return (int)_buf[_cur];
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until we see a newline, EOF, or until 'len' characters have been
	 * read.
	 */
	size_t gets(char *buf, size_t len) {
		size_t stored = 0;
		while(true) {
			int c = get();
			if(c == -1) {
				// End-of-file
				buf[stored] = '\0';
				return stored;
			}
			if(stored == len-1 || c == '\n' || c == '\r') {
				// End of string
				buf[stored] = '\0';
				// Skip over all end-of-line characters
				int pc = peek();
				while(pc == '\n' || pc == '\r') {
					get(); // discard
					pc = peek();
				}
				// Next get() will be after all newline characters
				return stored;
			}
			buf[stored++] = (char)c;
		}
	}

	/**
	 * Store a string of characters from the input file into 'buf',
	 * until we see a newline, EOF, or until 'len' characters have been
	 * read.
	 */
	size_t get(char *buf, size_t len) {
		size_t stored = 0;
		for(size_t i = 0; i < len; i++) {
			int c = get();
			if(c == -1) return i;
			buf[stored++] = (char)c;
		}
		return len;
	}

	static const size_t LASTN_BUF_SZ = 8 * 1024;

	/**
	 * Keep get()ing characters until a non-whitespace character (or
	 * -1) is reached, and return it.
	 */
	int getPastWhitespace() {
		int c;
		while(isspace(c = get()) && c != -1);
		return c;
	}

	/**
	 * Keep get()ing characters until a we've passed over the next
	 * string of newline characters (\r's and \n's) or -1 is reached,
	 * and return it.
	 */
	int getPastNewline() {
		int c = get();
		while(c != '\r' && c != '\n' && c != -1) c = get();
		while(c == '\r' || c == '\n') c = get();
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}

	/**
	 * Keep get()ing characters until a we've passed over the next
	 * string of newline characters (\r's and \n's) or -1 is reached,
	 * and return it.
	 */
	int peekPastNewline() {
		int c = peek();
		while(c != '\r' && c != '\n' && c != -1) c = get();
		while(c == '\r' || c == '\n') c = get();
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}

	/**
	 * Keep peek()ing then get()ing characters until the next return
	 * from peek() is just after the last newline of the line.
	 */
	int peekUptoNewline() {
		int c = peek();
		while(c != '\r' && c != '\n' && c != -1) {
			get(); c = peek();
		}
		while(c == '\r' || c == '\n') {
			get();
			c = peek();
		}
		assert_neq(c, '\r');
		assert_neq(c, '\n');
		return c;
	}

	size_t lastNCur() const { return _lastn_cur; }

	/**
	 * Reset to the beginning of the last-N-chars buffer.
	 */
	void resetLastN() {
		_lastn_cur = 0;
	}

	/**
	 * Copy the last several characters in the last-N-chars buffer
	 * (since the last reset) into the provided buffer.
	 */
	size_t copyLastN(char *buf) {
		memcpy(buf, _lastn_buf, _lastn_cur);
		return _lastn_cur;
	}

	/**
	 * Get const pointer to the last-N-chars buffer.
	 */
	const char *lastN() const {
		return _lastn_buf;
	}

	/**
	 * Get current size of the last-N-chars buffer.
	 */
	const size_t lastNLen() const {
		return _lastn_cur;
	}

private:

	void init() {
		_in = NULL;
		_inf = NULL;
		_ins = NULL;
		_cur = _buf_sz = BUF_SZ;
		_done = false;
		_lastn_cur = 0;
		// no need to clear _buf[]
	}

	static const size_t BUF_SZ = 256 * 1024;
	FILE     *_in;
	std::ifstream *_inf;
	std::istream  *_ins;
	size_t    _cur;
	size_t    _buf_sz;
	bool      _done;
	uint8_t   _buf[BUF_SZ]; // (large) input buffer
	size_t    _lastn_cur;
	char      _lastn_buf[LASTN_BUF_SZ]; // buffer of the last N chars dispensed
};

/**
 * Wrapper for a buffered output stream that writes bitpairs.
 */
class BitpairOutFileBuf {
public:
	/**
	 * Open a new output stream to a file with given name.
	 */
	BitpairOutFileBuf(const char *in) : bpPtr_(0), cur_(0) {
		assert(in != NULL);
		out_ = fopen(in, "wb");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open bitpair-output file " << in << std::endl;
			throw 1;
		}
		memset(buf_, 0, BUF_SZ);
	}

	/**
	 * Write a single bitpair into the buf.  Flush the buffer if it's
	 * full.
	 */
	void write(int bp) {
		assert_lt(bp, 4);
		assert_geq(bp, 0);
		buf_[cur_] |= (bp << bpPtr_);
		if(bpPtr_ == 6) {
			bpPtr_ = 0;
			cur_++;
			if(cur_ == BUF_SZ) {
				// Flush the buffer
				if(!fwrite((const void *)buf_, BUF_SZ, 1, out_)) {
					std::cerr << "Error writing to the reference index file (.4.ebwt)" << std::endl;
					throw 1;
				}
				// Reset to beginning of the buffer
				cur_ = 0;
			}
			// Initialize next octet to 0
			buf_[cur_] = 0;
		} else {
			bpPtr_ += 2;
		}
	}

	/**
	 * Write any remaining bitpairs and then close the input
	 */
	void close() {
		if(cur_ > 0 || bpPtr_ > 0) {
			if(bpPtr_ == 0) cur_--;
			if(!fwrite((const void *)buf_, cur_ + 1, 1, out_)) {
				std::cerr << "Error writing to the reference index file (.4.ebwt)" << std::endl;
				throw 1;
			}
		}
		fclose(out_);
	}
private:
	static const size_t BUF_SZ = 128 * 1024;
	FILE    *out_;
	int      bpPtr_;
	size_t   cur_;
	char     buf_[BUF_SZ]; // (large) input buffer
};

/**
 * Wrapper for a buffered output stream that writes characters and
 * other data types.  This class is *not* synchronized; the caller is
 * responsible for synchronization.
 */
class OutFileBuf {

public:

	/**
	 * Open a new output stream to a file with given name.
	 */
	OutFileBuf(const char *out, bool binary = false) :
		name_(out), cur_(0), closed_(false)
	{
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
		if(setvbuf(out_, NULL, _IOFBF, 10* 1024* 1024)) 
			std::cerr << "Warning: Could not allocate the proper buffer size for output file stream. " << std::endl;
	}

	/**
	 * Open a new output stream to standard out.
	 */
	OutFileBuf() : name_("cout"), cur_(0), closed_(false) {
		out_ = stdout;
	}

	/**
	 * Open a new output stream to a file with given name.
	 */
	void setFile(const char *out, bool binary = false) {
		assert(out != NULL);
		out_ = fopen(out, binary ? "wb" : "w");
		if(out_ == NULL) {
			std::cerr << "Error: Could not open alignment output file " << out << std::endl;
			throw 1;
		}
		reset();
	}

	/**
	 * Write a single character into the write buffer and, if
	 * necessary, flush.
	 */
	void write(char c) {
		assert(!closed_);
		if(cur_ == BUF_SZ) flush();
		buf_[cur_++] = c;
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	void writeString(const std::string& s) {
		assert(!closed_);
		size_t slen = s.length();
		if(cur_ + slen > BUF_SZ) {
			if(cur_ > 0) flush();
			if(slen >= BUF_SZ) {
				size_t wlen = fwrite(s.c_str(), 1, slen, out_);
				if(wlen != slen) {
					std::cerr << "Error while writing string output; " << slen
							  << " characters in string, " << wlen
							  << " written" << std::endl;
					throw 1;
				}
			} else {
				memcpy(&buf_[cur_], s.data(), slen);
				assert_eq(0, cur_);
				cur_ = slen;
			}
		} else {
			memcpy(&buf_[cur_], s.data(), slen);
			cur_ += slen;
		}
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write a c++ string to the write buffer and, if necessary, flush.
	 */
	void writeChars(const char * s, size_t len) {
		assert(!closed_);
		if(cur_ + len > BUF_SZ) {
			if(cur_ > 0) flush();
			if(len >= BUF_SZ) {
				size_t wlen = fwrite(s, 1, len, out_);
				if(wlen != len) {
					std::cerr << "Error while writing string output; " << len
							  << " characters in string, " << wlen
							  << " written" << std::endl;
					throw 1;
				}
			} else {
				memcpy(&buf_[cur_], s, len);
				assert_eq(0, cur_);
				cur_ = len;
			}
		} else {
			memcpy(&buf_[cur_], s, len);
			cur_ += len;
		}
		assert_leq(cur_, BUF_SZ);
	}

	/**
	 * Write any remaining bitpairs and then close the input
	 */
	void close() {
		if(closed_) return;
		if(cur_ > 0) flush();
		closed_ = true;
		if(out_ != stdout) {
			fclose(out_);
		}
	}

	/**
	 * Reset so that the next write is as though it's the first.
	 */
	void reset() {
		cur_ = 0;
		closed_ = false;
	}

	void flush() {
		if(cur_ != fwrite((const void *)buf_, 1, cur_, out_)){
			std::cerr << "Error while flushing and closing output" << std::endl;
			throw 1;
		}
		cur_ = 0;
	}

	/**
	 * Return true iff this stream is closed.
	 */
	bool closed() const {
		return closed_;
	}

	/**
	 * Return the filename.
	 */
	const char *name() {
		return name_;
	}

private:

	static const size_t BUF_SZ = 16 * 1024;

	const char *name_;
	FILE       *out_;
	size_t    cur_;
	char        buf_[BUF_SZ]; // (large) input buffer
	bool        closed_;
};

#endif /*ndef FILEBUF_H_*/
