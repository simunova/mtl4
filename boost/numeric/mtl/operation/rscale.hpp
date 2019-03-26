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
/*
 *  rscale.hpp
 *  MTL4
 *
 *  Created by Hui Li (huil@Princeton.EDU)
 *
 */

#ifndef MTL_RSCALE_INCLUDE
#define MTL_RSCALE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/utility/algebraic_category.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace tfunctor {
	
    // AlgebraicCategory is by default tag::scalar
    template <typename Value1, typename Value2, typename AlgebraicCategory>
    struct rscale
    {
	typedef typename Multiplicable<Value1, Value2>::result_type result_type;
	
	explicit rscale(const Value2& v2) : v2(v2) {}
	
	result_type operator() (const Value1& v1) const
	{
	    return v1 * v2;
	}

	Value2 value() const { return v2; }
    private:
	Value2 v2; 
    };
	
	
    template <typename Matrix, typename Value2>
    struct rscale<Matrix, Value2, tag::matrix>
    {
	typedef mat::rscaled_view<Matrix,Value2> result_type;
	
	explicit rscale(const Value2& v2) : v2(v2) {}
	
	result_type operator() (const Matrix& matrix) const
	{
	    return result_type(matrix, v2);
	}
    private:
	Value2 v2;
    };
	

    template <typename Vector, typename Value2>
    struct rscale<Vector, Value2, tag::vector>
    {
	typedef vec::rscaled_view<Vector, Value2> result_type;
	
	explicit rscale(const Value2& v2) : v2(v2) {}
		
	result_type operator() (const Vector& vector) const
	{
	    return result_type(vector, v2);
	}
    private:
	Value2 v2;
    };


} // namespace tfunctor
	

namespace mat {

    template <typename Value1, typename Value2>
    typename tfunctor::rscale<Value1, Value2, typename traits::algebraic_category<Value1>::type>::result_type
    inline rscale(const Value1& value1, const Value2& value2)
    {
	vampir_trace<4017> tracer;
	return tfunctor::rscale<Value1, Value2, typename traits::algebraic_category<Value1>::type>(value2)(value1);
    }
}

namespace vec {

    template <typename Value1, typename Value2>
    typename tfunctor::rscale<Value1, Value2, typename traits::algebraic_category<Value1>::type>::result_type
    inline rscale(const Value1& value1, const Value2& value2)
    {
	vampir_trace<2024> tracer;
	return tfunctor::rscale<Value1, Value2, typename traits::algebraic_category<Value1>::type>(value2)(value1);
    }
}

	
} // namespace mtl

#endif // MTL_RSCALE_INCLUDE
