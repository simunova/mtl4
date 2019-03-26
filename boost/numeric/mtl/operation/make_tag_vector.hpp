// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_MAKE_TAG_VECTOR_INCLUDE
#define MTL_MAKE_TAG_VECTOR_INCLUDE

#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/srange.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>

namespace mtl {

    template <typename Range>
    inline dense_vector<bool>
    make_tag_vector(std::size_t n, const Range& r)
    {
	dense_vector<bool> v(n);
	for (std::size_t i= 0; i < n; i++)
	    v[i]= r.in_range(i);
	return v;
    }

} // namespace mtl

#endif // MTL_MAKE_TAG_VECTOR_INCLUDE
