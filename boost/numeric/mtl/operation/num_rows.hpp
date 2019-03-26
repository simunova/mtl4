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

#ifndef MTL_NUM_ROWS_INCLUDE
#define MTL_NUM_ROWS_INCLUDE

#include <vector>

namespace mtl {

namespace traits {

    /// General declaration, used to disable unsupported types
    template <typename Collection>
    struct num_rows {};
    
    /// num_rows implementation for STL vectors
    template <typename Value>
    struct num_rows< std::vector<Value> > 
    {
	typedef std::size_t   type;
	type operator()(const std::vector<Value>& v) { return v.size(); }
    };

    /// num_rows implementation for (1D) arrays interpreted as vectors
    template <typename Value, unsigned Size>
    struct num_rows<Value[Size]>
    {
	typedef std::size_t   type;
	type operator()(const Value[Size]) { return Size; }
    };	   

    /// num_rows implementation for (2D and higher) arrays interpreted as matrices
    template <typename Value, unsigned Rows, unsigned Cols>
    struct num_rows<Value[Rows][Cols]>
    {
	typedef std::size_t   type;
	type operator()(const Value[Rows][Cols]) { return Rows; }
    };	    
}


/// num_rows function for non-MTL types (uses implicit enable_if), 1D interpreted as Column vector
template <typename Collection>
typename traits::num_rows<Collection>::type 
inline num_rows(const Collection& c)
{
    return traits::num_rows<Collection>()(c);
}


} // namespace mtl

#endif // MTL_NUM_ROWS_INCLUDE
