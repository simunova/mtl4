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

#ifndef MTL_SCALE_INCLUDE
#define MTL_SCALE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/true_copy.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { 

    namespace vec {
	
	template <typename Value1, typename Vector>
	typename traits::enable_if_vector<Vector, scaled_view<typename mtl::traits::true_copy<Value1>::type, Vector> >::type
	inline scale(const Value1& value1, const Vector& vector)
	{
	    vampir_trace<2028> tracer;
	    return scaled_view<typename mtl::traits::true_copy<Value1>::type, Vector>(value1, vector);
	}
    }

    namespace mat {

	template <typename Value1, typename Matrix>
	typename traits::enable_if_matrix<Matrix, scaled_view<Value1, Matrix> >::type
	inline scale(const Value1& value1, const Matrix& matrix)
	{
		vampir_trace<3029> tracer;
	    return scaled_view<Value1, Matrix> (value1, matrix);
	}
    }

    using vec::scale;
    using mat::scale;


    namespace tfunctor {

	// AlgebraicCategory is by default tag::scalar
	template <typename Value1, typename Value2, typename AlgebraicCategory>
	struct scale
	{
	    typedef typename Multiplicable<Value1, Value2>::result_type result_type;	    
	    explicit scale(const Value1& v1) : v1(v1) {}
	    
	    result_type operator() (const Value2& v2) const
	    {
		return v1 * v2;
	    }

	    Value1 value() const { return v1; }
	private:
	    Value1 v1; 
	};


	template <typename Value1, typename Matrix>
	struct scale<Value1, Matrix, tag::matrix>
	{
	    typedef mat::scaled_view<Value1, Matrix> result_type;	    
	    explicit scale(const Value1& v1) : v1(v1) {}
	
	    result_type operator() (const Matrix& matrix) const
	    {
		return result_type(v1, matrix);
	    }
	private:
	    typename mtl::traits::true_copy<Value1>::type v1; 
	};


	template <typename Value1, typename Vector>
	struct scale<Value1, Vector, tag::vector>
	{
	    typedef typename mtl::traits::true_copy<Value1>::type type1;

	    typedef vec::scaled_view<type1, Vector> result_type;
	    explicit scale(const Value1& v1) : v1(type1(v1)) {}
	
	    result_type operator() (const Vector& vector) const
	    {
		return result_type(v1, vector);
	    }
	private:
	    type1 v1; 
	};

    } // namespace tfunctor

} // namespace mtl

#endif // MTL_SCALE_INCLUDE
