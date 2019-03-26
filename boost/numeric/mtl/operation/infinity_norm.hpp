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

#ifndef MTL_INFINITY_NORM_INCLUDE
#define MTL_INFINITY_NORM_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/operation/max_of_sums.hpp>
#include <boost/numeric/mtl/vector/lazy_reduction.hpp>
#include <boost/numeric/mtl/vector/reduction.hpp>
#include <boost/numeric/mtl/vector/reduction_functors.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace vec {

	template <unsigned long Unroll, typename Vector>
	typename traits::enable_if_vector<Vector, typename RealMagnitude<typename Collection<Vector>::value_type>::type>::type
	inline infinity_norm(const Vector& vector)
	{
	    vampir_trace<2006> tracer;
	    typedef typename RealMagnitude<typename Collection<Vector>::value_type>::type result_type;
	    return reduction<Unroll, infinity_norm_functor, result_type>::apply(vector);
	}

	/*! Infinity-norm for vectors: infinity_norm(x) \f$\rightarrow |x|_\infty\f$.
	    \retval The magnitude type of the respective value type, see Magnitude.

	    The norms are defined as \f$|v|_\infty=\max_i |v_i|\f$.

	    Vector norms are unrolled 8-fold by default. 
	    An n-fold unrolling can be generated with infinity_norm<n>(x).
	    The maximum for n is 8 (it might be increased later).
	**/
	template <typename Vector>
	typename mtl::traits::enable_if_vector<Vector, typename RealMagnitude<typename Collection<Vector>::value_type>::type>::type
	inline infinity_norm(const Vector& vector)
	{
	    return infinity_norm<8>(vector);
	}

	template <typename Vector>
	lazy_reduction<Vector, infinity_norm_functor> inline lazy_infinity_norm(const Vector& v)
	{  return lazy_reduction<Vector, infinity_norm_functor>(v); 	}
    }

    namespace mat {
	
	// Ignore unrolling for matrices 
	template <unsigned long Unroll, typename Matrix>
	typename mtl::traits::enable_if_matrix<Matrix, typename RealMagnitude<typename Collection<Matrix>::value_type>::type>::type
	inline infinity_norm(const Matrix& matrix)
	{
	    vampir_trace<3011> tracer;
	    using mtl::impl::max_of_sums;
	    typename mtl::traits::row<Matrix>::type                             row(matrix); 
	    return max_of_sums(matrix, mtl::traits::is_row_major<typename OrientedCollection<Matrix>::orientation>(), 
			       row, num_rows(matrix));
	}

	/*! Infinity-norm for matrices: infinity_norm(x) \f$\rightarrow |x|_\infty\f$.
	    \retval The magnitude type of the respective value type, see Magnitude.

	    The norms are defined as \f$|A|_\infty=\max_i\{\sum_j(|A_{ij}|)\}\f$.
	    Matrix norms are not (yet) optimized by unrolling.
	**/
	template <typename Matrix>
	typename mtl::traits::enable_if_matrix<Matrix, typename RealMagnitude<typename Collection<Matrix>::value_type>::type>::type
	inline infinity_norm(const Matrix& matrix)
	{
	    return infinity_norm<8>(matrix);
	}
    }

    using vec::infinity_norm;
    using vec::lazy_infinity_norm;
    using mat::infinity_norm;

} // namespace mtl

#endif // MTL_INFINITY_NORM_INCLUDE
