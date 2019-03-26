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

#ifndef MTL_MAX_INCLUDE
#define MTL_MAX_INCLUDE

#include <iostream>
#include <cmath>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/vector/reduction.hpp>
#include <boost/numeric/mtl/vector/reduction_functors.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace vec {

    namespace impl {

	// Do we really need this for matrices?
	// Then we need a different dispatching
	
	template <unsigned long Unroll, typename Vector>
	typename Collection<Vector>::value_type
	inline max(const Vector& vector, tag::vector)
	{
	    typedef typename Collection<Vector>::value_type result_type;
	    return vec::reduction<Unroll, vec::max_functor, result_type>::apply(vector);
	}
	
    } // namespace impl

///Returns maximal entry of %vector v
template <unsigned long Unroll, typename Value>
typename Collection<Value>::value_type
inline max(const Value& value)
{
    vampir_trace<2010> tracer;
    return impl::max<Unroll>(value, typename traits::category<Value>::type());
}

template <typename Value>
typename Collection<Value>::value_type
inline max(const Value& value)
{
    return max<8>(value);
}

} // namespace vector

using vec::max;

} // namespace mtl::vector

#endif // MTL_MAX_INCLUDE
