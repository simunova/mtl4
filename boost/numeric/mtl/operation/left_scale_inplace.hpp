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

#ifndef MTL_LEFT_SCALE_INPLACE_INCLUDE
#define MTL_LEFT_SCALE_INPLACE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/operation/assign_each_nonzero.hpp>
#include <boost/numeric/mtl/operation/mult.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace impl {

	template <typename Factor, typename Coll>
	inline Coll& left_scale_inplace(const Factor& alpha, tag::scalar, Coll& c)
	{
	    vampir_trace<3014> tracer;	    
	    assign_each_nonzero(c, tfunctor::scale<Factor, typename Collection<Coll>::value_type>(alpha));

	    return c;
	}

	template <typename Factor, typename Collection>
	inline Collection& left_scale_inplace(const Factor& alpha, tag::matrix, Collection& c)
	{
	    vampir_trace<3014> tracer;	    
	    Collection tmp(alpha * c);
	    swap(tmp, c);
	    return c;
	}
    }
    
    namespace mat {
	
	/// Scale matrix \p c from left with scalar or matrix factor \p alpha; \p c is altered
	template <typename Factor, typename Matrix>
	typename mtl::traits::enable_if_matrix<Matrix, Matrix&>::type
	inline left_scale_inplace(const Factor& alpha, Matrix& A)
	{
	    return mtl::impl::left_scale_inplace(alpha, typename mtl::traits::category<Factor>::type(), A);
	}
    }

    namespace vec {
	
	/// Scale vector \p c from left with scalar or matrix factor \p alpha; \p c is altered
	template <typename Factor, typename Vector>
	typename mtl::traits::enable_if_vector<Vector, Vector&>::type
	inline left_scale_inplace(const Factor& alpha, Vector& v)
	{
	    return mtl::impl::left_scale_inplace(alpha, typename mtl::traits::category<Factor>::type(), v);
	}
    }

    using vec::left_scale_inplace;
    using mat::left_scale_inplace;

} // namespace mtl

#endif // MTL_LEFT_SCALE_INPLACE_INCLUDE
