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

#ifndef MTL_SIZE_INCLUDE
#define MTL_SIZE_INCLUDE

#include <vector>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

namespace traits {

    /// General declaration, used to disable unsupported types
    template <typename Collection>
    struct size {};
    
    /// size implementation for STL vectors
    template <typename Value>
    struct size< std::vector<Value> > 
    {
	typedef std::size_t   type;
	type operator()(const std::vector<Value>& v) { return v.size(); }
    };

    /// size implementation for (1D) arrays interpreted as vectors
    template <typename Value, unsigned Size>
    struct size<Value[Size]>
    {
    vampir_trace<2031> tracer;
	typedef std::size_t   type;
	type operator()(const Value[Size]) { return Size; }
    };	   

    /// size implementation for (2D and higher) arrays interpreted as matrices
    template <typename Value, unsigned Rows, unsigned Cols>
    struct size<Value[Rows][Cols]>
    {
	typedef std::size_t   type;
	vampir_trace<3033> tracer;
	type operator()(const Value[Rows][Cols]) { return Rows * Cols; }
    };	    
}


/// size function for non-MTL types (uses implicit enable_if)
template <typename Collection>
typename traits::size<Collection>::type 
inline size(const Collection& c)
{
	
    return traits::size<Collection>()(c);
}


} // namespace mtl

#endif // MTL_SIZE_INCLUDE
