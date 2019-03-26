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

#ifndef MTL_RESOURCE_INCLUDE
#define MTL_RESOURCE_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace traits {

	template <typename Vector>
	struct vector_resource
	{
	    typedef typename Collection<Vector>::size_type type;
	    type inline static apply(const Vector& v) { using mtl::vec::size; return mtl::size(v); }
	};
    }

    namespace vec {

	/// Describes the resources need for a certain vector.
	/** All necessary information to construct appropriate/consistent temporary vectors.
	    Normally, this is just the size of the vector.
	    For distributed vector we also need its distribution. **/
	template <typename Vector>
	typename mtl::traits::vector_resource<Vector>::type
	inline resource(const Vector& v)
	{
	    vampir_trace<4> tracer;
	    return mtl::traits::vector_resource<Vector>::apply(v);
	}

    } // namespace vector

    namespace mat {
	// maybe a pair of size_type? like position
    }

} // namespace mtl

#endif // MTL_RESOURCE_INCLUDE
