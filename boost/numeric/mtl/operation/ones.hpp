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

#ifndef MTL_ONES_INCLUDE
#define MTL_ONES_INCLUDE

#include <boost/numeric/mtl/matrix/implicit_dense.hpp>

namespace mtl {

    namespace mat {

	/// Return r by c matrix with ones of type Value in all entries
	template <typename Value>
	ones_matrix<Value> inline ones(std::size_t r, std::size_t c)
	{
	    return ones_matrix<Value>(r, c);
	}

	/// Return r by c matrix with ones of type Value in all entries
	ones_matrix<> inline ones(std::size_t r, std::size_t c)
	{
	    return ones_matrix<>(r, c);
	}

    }

    using mat::ones;


} // namespace mtl

#endif // MTL_ONES_INCLUDE
