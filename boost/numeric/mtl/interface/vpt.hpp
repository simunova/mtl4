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

#ifndef MTL_VPT_VPT_INCLUDE
#define MTL_VPT_VPT_INCLUDE

#ifdef MTL_HAS_VPT
  #include <vt_user.h> 
  #include <boost/mpl/bool.hpp>
#endif 

#include <math.h> 
#include <string>

namespace mtl { 

/// Namespace for Vampir Trace interface
namespace vpt {

#ifdef MTL_HAS_VPT

#ifndef MTL_VPT_LEVEL
#  define MTL_VPT_LEVEL 2
#endif 

/// Class for Vampir Trace
template <int N>
class vampir_trace
{
    // Statically determine whether the event is traced; just in case you wanted to know how.
    typedef boost::mpl::bool_<(MTL_VPT_LEVEL * 1000 < N)> to_print;
  public:
    /// Default constructor defines the start point of a trace
    vampir_trace() { entry(to_print());  }

    void entry(boost::mpl::false_) {}
    void entry(boost::mpl::true_) 
    {	
	VT_USER_START(name.c_str()); 
	// std::cout << "vpt_entry(" << N << ")\n";    
    }
    
    /// Destructor defines the end point of a trace
    ~vampir_trace() { end(to_print());  }

    void end(boost::mpl::false_) {}
    void end(boost::mpl::true_) 
    {
	VT_USER_END(name.c_str()); 
	// std::cout << "vpt_end(" << N << ")\n";    
    }
    
    /// Function to check whether this event is traced with the current setting
    bool is_traced() { return to_print::value; }

  private:
    static std::string name;
};


#else

// Dummy when Vampir Trace is not supported
template <int N>
class vampir_trace 
{
  public:
    vampir_trace() {}
    void show_vpt_level() {}
    bool is_traced() { return false; }
  private:
    static std::string name;
};
#endif

    // names defined in vpt.cpp !!!

} // namespace vpt

/// Import of vpt::vampir_trace
using vpt::vampir_trace;

} // namespace mtl

#endif // MTL_VPT_VPT_INCLUDE
