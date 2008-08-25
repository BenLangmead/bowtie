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
  $Id: system_sema.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_SEMAPHORE_H
#define SEQAN_HEADER_SYSTEM_SEMAPHORE_H

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES SemaphoreDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Semaphore
    {
        typedef LONG Type;
        typedef HANDLE Handle;
        enum { MAX_VALUE = MAXLONG };
        
        Handle hSemaphore;

        Semaphore(Type init = 0, Type max = MAX_VALUE) {
            SEQAN_DO_SYS2((hSemaphore = CreateSemaphore(&SemaphoreDefaultAttributes, init, max, NULL)) != NULL, "Could not create Semaphore")
        }

        ~Semaphore() {
            SEQAN_DO_SYS2(CloseHandle(hSemaphore) != 0, "Could not destroy Semaphore")
        }

        bool lock(DWORD timeout_millis = INFINITE) {
            return WaitForSingleObject(hSemaphore, timeout_millis) != WAIT_TIMEOUT;
        }

        void unlock() {
            SEQAN_DO_SYS2(ReleaseSemaphore(hSemaphore, 1, NULL) != 0, "Could not unlock Semaphore")
        }

    private:

        Semaphore(Semaphore const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };

#else

    struct Semaphore
    {
        typedef unsigned int Type;
        typedef sem_t* Handle;
        
        sem_t data, *hSemaphore;

        Semaphore(Type init = 0):
            hSemaphore(&data)
        {
            SEQAN_DO_SYS(!sem_init(hSemaphore, 0, init));
        }

        ~Semaphore() {
            SEQAN_DO_SYS(!sem_destroy(hSemaphore));
        }

        void lock() {
            SEQAN_DO_SYS(!sem_wait(hSemaphore));
        }

        void unlock() {
            SEQAN_DO_SYS(!sem_post(hSemaphore));
        }

    private:

        Semaphore(Semaphore const &) {
            // resource copying is not yet supported (performance reason)
            // it needs a reference counting technique
        }

    };


#endif

}

#endif
