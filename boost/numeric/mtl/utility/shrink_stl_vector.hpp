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

#ifndef MTL_SHRINK_STL_VECTOR_INCLUDE
#define MTL_SHRINK_STL_VECTOR_INCLUDE

namespace mtl {

/// Shrink memory consumption of an STL vector to its size
template <typename Value, typename Allocator>
void inline shrink_stl_vector(std::vector<Value, Allocator>& v)
{
    if (v.capacity() > v.size()) {
	std::vector<Value, Allocator> tmp(v.begin(), v.end());
	swap(tmp, v);
    }
}


} // namespace mtl

#endif // MTL_SHRINK_STL_VECTOR_INCLUDE
