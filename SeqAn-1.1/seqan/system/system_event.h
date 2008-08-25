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
  $Id: system_event.h,v 1.1 2008/08/25 16:20:08 langmead Exp $
 ==========================================================================*/

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SYSTEM_EVENT_H
#define SEQAN_HEADER_SYSTEM_EVENT_H

//#include <iterator>

namespace SEQAN_NAMESPACE_MAIN
{

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES EventDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

    struct Event    // this class mustn't exceed the size of HANDLE (needed by waitForAll/Any)
    {
        typedef HANDLE Handle;
		enum { Infinite = INFINITE };
        Handle hEvent;

        Event():
            hEvent(NULL) {}

        Event(BOOL initial) {
            SEQAN_DO_SYS2(open(initial), "Could not create Event")
        }

        ~Event() {
            if (*this) SEQAN_DO_SYS2(close(), "Could not destroy Event")
        }

        Event(Event const &origin) {
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                hEvent = origin.hEvent;
                const_cast<Event&>(origin).hEvent = NULL;
            } else
                hEvent = NULL;
        }

        inline Event& operator=(Event const &origin) {
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                hEvent = origin.hEvent;
                const_cast<Event&>(origin).hEvent = NULL;
            } else
                hEvent = NULL;
            return *this;
        }

        inline bool open(BOOL initial = FALSE) {
            return (hEvent = CreateEvent(&EventDefaultAttributes, TRUE, initial, NULL)) != NULL;
        }

        inline bool close() {
            return CloseHandle(hEvent) && !(hEvent = NULL);
        }

        inline bool wait(DWORD timeout_millis = Infinite) {
            if (!hEvent) return true;
            return WaitForSingleObject(hEvent, timeout_millis) != WAIT_TIMEOUT;
        }

        inline bool signal() {
            return SetEvent(hEvent) != 0;
        }

        inline bool reset() {
            return ResetEvent(hEvent) != 0;
        }

        inline operator bool() const {
            return hEvent != NULL;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global event functions
	
	inline void reset(Event &e) {
		e.reset();
	}

    inline bool waitForAll(Event eventList[], DWORD count, DWORD timeout_millis)
	{
		return WaitForMultipleObjects(count, &eventList[0].hEvent, true, timeout_millis) != WAIT_TIMEOUT;
	}
    
    inline bool waitForAll(Event eventList[], DWORD count)
	{
		return WaitForMultipleObjects(count, &eventList[0].hEvent, true, Event::Infinite) != WAIT_TIMEOUT;
	}
    
	inline int waitForAny(Event eventList[], DWORD count, DWORD timeout_millis)
	{
        DWORD result = WaitForMultipleObjects(count, &eventList[0].hEvent, false, timeout_millis = Event::Infinite);

        if (/*result >= WAIT_OBJECT_0 && */result < WAIT_OBJECT_0 + count)
    		return result - WAIT_OBJECT_0;

        return -1;
	}
    
#else

    struct Event: public Mutex
    {
        typedef pthread_cond_t* Handle;
		enum { Infinite = LONG_MAX };
        pthread_cond_t data, *hEvent;

        Event():
            hEvent(NULL) {}

        Event(bool initial) {
            SEQAN_DO_SYS(open(initial));
        }

        ~Event() {
            if (*this)
                SEQAN_DO_SYS(close());
        }

		Event(Event const &origin):
			Mutex()
		{
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                data = origin.data;
                const_cast<Event&>(origin).hEvent = NULL;
                hEvent = &data;
            } else
                hEvent = NULL;
        }

        inline Event& operator=(Event const &origin) {
            // resource sharing is not yet supported (performance reason)
            // it needs a reference counting technique
            if (origin) {
                data = origin.data;
                const_cast<Event&>(origin).hEvent = NULL;
                hEvent = &data;
            } else
                hEvent = NULL;
            return *this;
        }

        inline bool open(bool initial = false)
        {
            if (Mutex::open() && !pthread_cond_init(&data, NULL) && (hEvent = &data)) {
                if (initial) return signal();
                return true;
            } else
                return false;
        }

        inline bool close() {
            return !(pthread_cond_destroy(hEvent) || (hEvent = NULL) || !Mutex::close());
        }

        inline bool wait() {
            if (!hEvent) return true;
            Mutex::lock();
            return !pthread_cond_wait(hEvent, Mutex::hMutex);
        }

        inline bool wait(long timeout_millis) {
            if (timeout_millis != Infinite) {
                timespec ts;
                ts.tv_sec = timeout_millis / 1000;
                ts.tv_nsec = (timeout_millis % 1000) * 1000;
                Mutex::lock();
				return pthread_cond_timedwait(hEvent, Mutex::hMutex, &ts) != ETIMEDOUT;
            } else
                return wait();
        }

        inline bool signal() {
            return !pthread_cond_broadcast(hEvent);
        }

        inline operator bool() const {
            return hEvent != NULL;
        }
    };
    
#endif


	//////////////////////////////////////////////////////////////////////////////
	// global event functions

	inline bool open(Event &e, bool initial) {
		return e.open(initial);
	}

	inline bool open(Event &e) {
		return open(e, false);
	}

	inline bool close(Event &e) {
		return e.close();
	}

	inline bool waitFor(Event &e) {
		return e.wait();
	}

    template < typename TTime >
	inline bool waitFor(Event &e, TTime timeout_millis) {
        #ifdef SEQAN_PROFILE
            double begin = sysTime();
		    bool b = e.wait(timeout_millis);
            double end = sysTime();
            if (begin != end)
                ::std::cerr << "waitTime: " << end - begin << ::std::endl;
            return b;
        #else
            return e.wait(timeout_millis);
        #endif
	}

	inline bool signal(Event &e) {
		return e.signal();
	}

/*
	//////////////////////////////////////////////////////////////////////////////
	// emulate events in a singlethreaded environment

	struct DummyEvent {
		typedef	void Handle;
		DummyEvent(bool initial = false) {}
		inline bool wait(unsigned timeOut = NULL) { return true; }
		inline void reset() {}
		inline void signal() {}
	};

	//////////////////////////////////////////////////////////////////////////////
	// global dummy event functions

	template < typename TCount >
	inline bool waitForAll(DummyEvent eventList[], TCount count) {
		return true;
	}

	template < typename TCount, typename TTime >
	inline bool waitForAll(DummyEvent eventList[], TCount count, TTime timeout_millis) {
		return true;
	}

	template < typename TCount >
	inline TCount waitForAny(DummyEvent eventList[], TCount count) {
		return 0;
	}
	template < typename TCount, typename TTime >
	inline TCount waitForAny(DummyEvent eventList[], TCount count, TTime timeout_millis) {
		return 0;
	}
*/
}

#endif
