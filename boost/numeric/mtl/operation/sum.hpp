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

#ifndef MTL_SUM_INCLUDE
#define MTL_SUM_INCLUDE

#include <iostream>
#include <cmath>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/vector/lazy_reduction.hpp>
#include <boost/numeric/mtl/vector/reduction.hpp>
#include <boost/numeric/mtl/vector/reduction_functors.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace impl {

	// Do we really need this for matrices?
	
	template <unsigned long Unroll, typename Vector>
	typename Collection<Vector>::value_type
	inline sum(const Vector& vector, tag::vector)
	{
		vampir_trace<2035> tracer;
	    typedef typename Collection<Vector>::value_type result_type;
	    return vec::reduction<Unroll, vec::sum_functor, result_type>::apply(vector);
	}
	
    } // namespace impl

///Return sum of all %vector-entries
template <unsigned long Unroll, typename Value>
typename Collection<Value>::value_type
inline sum(const Value& value)
{
    return impl::sum<Unroll>(value, typename traits::category<Value>::type());
}

template <typename Value>
typename Collection<Value>::value_type
inline sum(const Value& value)
{
    return sum<8>(value);
}

namespace vec {
	template <typename Vector>
	lazy_reduction<Vector, sum_functor> inline lazy_sum(const Vector& v)
	{  return lazy_reduction<Vector, sum_functor>(v); 	}
}

using vec::lazy_sum;

} // namespace mtl

#endif // MTL_SUM_INCLUDE
