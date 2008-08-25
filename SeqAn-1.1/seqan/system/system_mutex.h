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
  $Id: system_mutex.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_MUTEX_H
#define SEQAN_HEADER_SYSTEM_MUTEX_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES MutexDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Mutex
    {
        typedef HANDLE Handle;

        Handle hMutex;

        Mutex():
            hMutex(NULL) {}

        Mutex(BOOL initial) {
            SEQAN_DO_SYS2(open(initial), "Could not create Mutex")
        }

        ~Mutex() {
            if (*this)
                SEQAN_DO_SYS2(close(), "Could not destroy Mutex")
        }

        inline bool open(BOOL initial = false) {
            return (hMutex = CreateMutex(&MutexDefaultAttributes, initial, NULL)) != NULL;
        }

        inline bool close() {
            return CloseHandle(hMutex) && !(hMutex = NULL);
        }

        inline bool lock(DWORD timeout_millis = INFINITE) {
            return WaitForSingleObject(hMutex, timeout_millis) != WAIT_TIMEOUT;
        }

        inline bool unlock() {
            return ReleaseMutex(hMutex) != 0;
        }

        inline operator bool() const {
            return hMutex != NULL;
        }

    private:

        Mutex(Mutex const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }
    };
    
#else

    struct Mutex
    {
        typedef pthread_mutex_t* Handle;
        
        pthread_mutex_t data, *hMutex;

        Mutex():
            hMutex(NULL) {}

        Mutex(bool initial) {
            SEQAN_DO_SYS(open(initial));
        }

        ~Mutex() {
            if (*this)
                SEQAN_DO_SYS(close());
        }

        inline bool open(bool initial = false)
        {
            if (!pthread_mutex_init(&data, NULL) && (hMutex = &data)) {
                if (initial) return lock();
                return true;
            } else
                return false;
        }

        inline bool close() {
            return !(pthread_mutex_destroy(hMutex) || (hMutex = NULL));
        }

        inline bool lock() {
            return !pthread_mutex_lock(hMutex);
        }

        inline bool unlock() {
            return !pthread_mutex_unlock(hMutex);
        }

        inline operator bool() const {
            return hMutex != NULL;
        }

    private:

        Mutex(Mutex const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };
    
#endif


	//////////////////////////////////////////////////////////////////////////////
	// global mutex functions

	inline bool open(Mutex &m, bool initial) {
		return m.open(initial);
	}

	inline bool open(Mutex &m) {
		return open(m, false);
	}

	inline bool close(Mutex &m) {
		return m.close();
	}

	inline bool lock(Mutex &m) {
		return m.lock();
	}

	inline bool unlock(Mutex &m) {
		return m.unlock();
	}

}

#endif
