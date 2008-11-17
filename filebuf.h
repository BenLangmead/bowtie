#ifndef FILEBUF_H_
#define FILEBUF_H_

#include <iostream>
#include <stdio.h>
#include <stdint.h>

/**
 * Simple wrapper for a FILE* that reads it in chunks (with fread) and
 * keeps those chunks in a buffer.  It also services calls to get(),
 * peek() and gets() from the buffer, reading in additional chunks when
 * necessary.  TODO: Implement asynchronous I/O.  Obstacle: neither
 * Cygwin nor MinGW seem to support POSIX aio.
 */
class FileBuf {
public:
	FileBuf(FILE *__in) : _in(__in), _ins(NULL), _cur(BUF_SZ), _buf_sz(BUF_SZ), _done(false) {
		assert(_in != NULL);
	}

	FileBuf(istream *__ins) : _in(NULL), _ins(__ins), _cur(BUF_SZ), _buf_sz(BUF_SZ), _done(false) {
		assert(_ins != NULL);
	}

	int get() {
		assert(_in != NULL || _ins != NULL);
		int c = peek();
		if(c != -1) _cur++;
		return c;
	}

	bool eof() {
		return (_cur == _buf_sz) && _done;
	}

	void reset() {
		if(_ins != NULL) {
			_ins->clear();
			_ins.seekg(0, ios::beg);
		} else {
			rewind(_in);
		}
		_cur = BUF_SZ;
		_buf_sz = BUF_SZ;
		_done = false;
	}

	int peek() {
		assert(_in != NULL || _ins != NULL);
		assert_leq(_cur, _buf_sz);
		if(_cur == _buf_sz) {
			if(_done) { return -1; }
			else {
				// Get the next chunk
				if(_ins != NULL) {
					_ins->read(_buf, BUF_SZ);
					_buf_sz = _ins->gcount();
				} else {
					assert(_in != NULL);
					_buf_sz = fread(_buf, 1, BUF_SZ, _in);
				}
				_cur = 0;
				if(_buf_sz == 0) {
					_done = true;
					return -1;
				} else if(_buf_sz < BUF_SZ) {
					_done = true;
				}
			}
		}
		return (int)_buf[_cur];
	}

	size_t gets(char *buf, size_t len) {
		size_t stored = 0;
		while(true) {
			int c = get();
			if(c == -1) {
				buf[stored] = '\0';
				return stored;
			}
			if(stored == len-1 || c == '\n' || c == '\r') {
				buf[stored] = '\0';
				// Skip over all end-of-line characters
				int pc = peek();
				while(pc == '\n' || pc == '\r') {
					get();
					pc = peek();
				}
				// Next get() will be after all newline characters
				return stored;
			}
			buf[stored++] = (char)c;
		}
	}
private:
	static const size_t BUF_SZ = 256 * 1024;
	FILE    *_in;
	istream *_ins;
	size_t   _cur;
	size_t   _buf_sz;
	bool     _done;
	char     _buf[BUF_SZ]; // (large) input buffer
};

#endif /*ndef FILEBUF_H_*/
