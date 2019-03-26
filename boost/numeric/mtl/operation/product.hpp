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

#ifndef MTL_PRODUCT_INCLUDE
#define MTL_PRODUCT_INCLUDE

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

    namespace vec {
	
	namespace impl {

	    // Do we really need this for matrices?
	
	    template <unsigned long Unroll, typename Vector>
	    typename Collection<Vector>::value_type
	    inline product(const Vector& vector, tag::vector)
	    {
		typedef typename Collection<Vector>::value_type result_type;
		return vec::reduction<Unroll, vec::product_functor, result_type>::apply(vector);
	    }
	
	} // namespace impl

	///Returns product of all collection-entries (%vector-entries)
	template <unsigned long Unroll, typename Value>
	typename Collection<Value>::value_type
	inline product(const Value& value)
	{	vampir_trace<2020> tracer;
	    return impl::product<Unroll>(value, typename traits::category<Value>::type());
	}

	template <typename Value>
	typename Collection<Value>::value_type
	inline product(const Value& value)
	{
	    return product<8>(value);
	}

	template <typename Vector>
	lazy_reduction<Vector, product_functor> inline lazy_product(const Vector& v)
	{  return lazy_reduction<Vector, product_functor>(v); 	}

    } // namespace vector

    using vec::lazy_product;
    using vec::product;

} // namespace mtl

#endif // MTL_PRODUCT_INCLUDE
