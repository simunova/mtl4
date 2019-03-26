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

#ifndef MTL_CROP_INCLUDE
#define MTL_CROP_INCLUDE

#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

    namespace vec {

	/// Remove all zero entries from a collection
	/** Does nothing for dense collections **/
	template <typename T>
	typename mtl::traits::enable_if_vector<T, T&>::type inline crop(T& x)
	{
	    vampir_trace<3006> tracer;
	    x.crop(); return x;
	}
    }

    namespace mat {

	/// Remove all zero entries from a collection
	/** Does nothing for dense collections **/
	template <typename T>
	typename mtl::traits::enable_if_matrix<T, T&>::type inline crop(T& x)
	{
	    vampir_trace<3006> tracer;
	    x.crop(); return x;
	}
    }
    
    using vec::crop;
    using mat::crop;

} // namespace mtl

#endif // MTL_CROP_INCLUDE
