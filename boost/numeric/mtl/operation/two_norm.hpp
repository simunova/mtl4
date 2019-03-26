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

#ifndef MTL_TWO_NORM_INCLUDE
#define MTL_TWO_NORM_INCLUDE

#include <iostream>
#include <cmath>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/vector/lazy_reduction.hpp>
#include <boost/numeric/mtl/vector/reduction.hpp>
#include <boost/numeric/mtl/vector/reduction_functors.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace vec {

	template <unsigned long Unroll, typename Value>
	typename RealMagnitude<typename Collection<Value>::value_type>::type
	inline two_norm(const Value& value)
	{
	    using std::sqrt;
	    vampir_trace<2039> tracer;
	    typedef typename RealMagnitude<typename Collection<Value>::value_type>::type result_type;
		return sqrt(reduction<Unroll, two_norm_functor, result_type>::apply(value));
	}

	
	/*! Two-norm for vectors: two_norm(x) \f$\rightarrow |x|_2\f$.
	    \retval The magnitude type of the respective value type, see Magnitude.
	    The norms are defined as \f$|v|_2=\sqrt{\sum_i |v_i|^2}\f$.
	    
	    Vector norms are unrolled 8-fold by default. 
	    An n-fold unrolling can be generated with two_norm<n>(x).
	    The maximum for n is 8 (it might be increased later).
	**/
	template <typename Value>
	typename RealMagnitude<typename Collection<Value>::value_type>::type
	inline two_norm(const Value& value)
	{
	    return two_norm<4>(value);
	}

	template <typename Vector>
	lazy_reduction<Vector, two_norm_functor> inline lazy_two_norm(const Vector& v)
	{  return lazy_reduction<Vector, two_norm_functor>(v); 	}
	

    } // namespace vector

    // two_norm for matrices not implemented (would need enable_if like one_norm)

    using vec::two_norm;
    using vec::lazy_two_norm;

} // namespace mtl

#endif // MTL_TWO_NORM_INCLUDE
