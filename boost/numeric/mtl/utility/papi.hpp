// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef MTL_PAPI_INCLUDE
#define MTL_PAPI_INCLUDE

#ifdef MTL_HAS_PAPI

#include <papi.h>
#include <boost/numeric/mtl/utility/exception.hpp>

#endif // MTL_HAS_PAPI

namespace mtl { namespace utility {

/// Exception for errors with PAPI, is sub-divided further
struct papi_error : public runtime_error
{
    explicit papi_error(const char *s= "PAPI error") : runtime_error(s) {}
};

struct papi_version_mismatch : public papi_error 
{
    explicit papi_version_mismatch(const char *s= "PAPI: version mismatch") : papi_error(s) {}
};

struct papi_no_counters : public papi_error 
{
    explicit papi_no_counters(const char *s= "PAPI: no counters") : papi_error(s) {}
};

struct papi_create_eventset_error : public papi_error 
{
    explicit papi_create_eventset_error(const char *s= "PAPI: create event set error") : papi_error(s) {}
};

struct papi_name_to_code_error : public papi_error 
{
    explicit papi_name_to_code_error(const char *s= "PAPI: name to code error") : papi_error(s) {}
};

struct papi_query_event_error : public papi_error 
{
    explicit papi_query_event_error(const char *s= "PAPI: query event error") : papi_error(s) {}
};

struct papi_start_event_error : public papi_error
{
    explicit papi_start_event_error(const char *s= "PAPI: start event error") : papi_error(s) {}
};

struct papi_add_event_error : public papi_error 
{
    explicit papi_add_event_error(const char *s= "PAPI: add event error") : papi_error(s) {}
};

struct papi_reset_error : public papi_error 
{
    explicit papi_reset_error(const char *s= "PAPI: reset error") : papi_error(s) {}
};

struct papi_read_error : public papi_error 
{
    explicit papi_read_error(const char *s= "PAPI: read error") : papi_error(s) {}
};

struct papi_index_range_error : public papi_error 
{
    explicit papi_index_range_error(const char *s= "PAPI: index range error") : papi_error(s) {}
};


#ifdef MTL_HAS_PAPI

class papi_t 
{
    void init_papi()
    {
	static bool initialized= false;
	if (!initialized) {

	    MTL_THROW_IF(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT,
			 papi_version_mismatch());

	    num_counters = PAPI_get_opt(PAPI_MAX_HWCTRS, NULL);
	    MTL_THROW_IF(num_counters <= 0, papi_no_counters());

	    counters= new long_long[num_counters];

	    MTL_THROW_IF(PAPI_create_eventset(&event_set) != PAPI_OK, papi_create_eventset_error());
	    initialized= true;
	}
    }

public:
    const static bool true_papi = true;

    papi_t() : event_set(PAPI_NULL), active_events(0)
    {
	init_papi();
    }


    ~papi_t()
    {
	delete[](counters);
    }


    // returns index of added event
    int add_event(const char* name)
    {
	int code;
	MTL_THROW_IF(PAPI_event_name_to_code(const_cast<char*>(name), &code) != PAPI_OK, 
		     papi_name_to_code_error());
	// std::cout << "add event " << const_cast<char*>(name) << " " << code << "\n";
	MTL_THROW_IF (PAPI_query_event(code) != PAPI_OK,
		      papi_query_event_error());
	MTL_THROW_IF (PAPI_add_event(event_set, code) != PAPI_OK,
		      papi_add_event_error());
	list_events();
	return active_events++;
    }

    void start() 
    {
	MTL_THROW_IF (PAPI_start(event_set) != PAPI_OK,
		      papi_start_event_error());
	reset();
    }

    void list_events()
    {
#if 0
	int evv[8], num= 8;
	PAPI_list_events(event_set, evv, &num);
	for (int i= 0; i < num; i++ ) std::cout << evv[i] << "\n";
#endif
    }

    bool is_event_supported(const char* name) const
    {
	int code;
	return PAPI_event_name_to_code(const_cast<char*>(name), &code) == PAPI_OK 
	       && PAPI_query_event(code) == PAPI_OK;
    }

    void reset()
    {
	MTL_THROW_IF(PAPI_reset(event_set) != PAPI_OK, papi_reset_error());
    }

    void read()
    {
	list_events();
	MTL_THROW_IF (PAPI_read(event_set, counters) != PAPI_OK, papi_read_error());
	// std::cout << "counters read, first value: " << *counters << "\n";
    }

    long_long operator[](int index) const
    {
	if (index < 0 || index >= active_events) throw papi_index_range_error();
	return counters[index];
    }

    //private:
    int num_counters, event_set, active_events;
    long_long *counters;
};

#else // no papi

// Faked papi type:

struct papi_t
{
    const static bool true_papi = false;
    int add_event(const char* name) { return 0;}
    void start() {}
    bool is_event_supported(const char* name) const { return false;}
    void reset() {}
    void read() {}
    long long operator[](int index) const { return 0; }
};
 
#endif

}} // namespace mtl::utility


#endif // MTL_PAPI_INCLUDE
