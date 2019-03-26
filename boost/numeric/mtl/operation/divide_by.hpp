/*
 *  divide_by.h
 *  MTL4
 *
 *  Created by Hui Li (huil@Princeton.EDU)
 *
 */

#ifndef MTL_DIVIDE_BY_INCLUDE
#define MTL_DIVIDE_BY_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/utility/algebraic_category.hpp>

namespace mtl { namespace tfunctor {
	
    // AlgebraicCategory is by default tag::scalar
    template <typename Value1, typename Value2, typename AlgebraicCategory>
    struct divide_by
    {
	typedef typename Multiplicable<Value1, Value2>::result_type result_type;
	
	explicit divide_by(const Value2& v2) : v2(v2) {}
	
	result_type operator() (const Value1& v1) const
	{
	    return v1 / v2;
	}
    private:
	Value2 v2; 
    };
	
	
    template <typename Matrix, typename Value2>
    struct divide_by<Matrix, Value2, tag::matrix>
    {
	typedef mat::divide_by_view<Matrix,Value2> result_type;
	
	explicit divide_by(const Value2& v2) : v2(v2) {}
	
	result_type operator() (const Matrix& matrix) const
	{
	    return result_type(matrix, v2);
	}
    private:
	Value2 v2;
    };
	

    template <typename Vector, typename Value2>
    struct divide_by<Vector, Value2, tag::vector>
    {
	typedef vec::divide_by_view<Vector, Value2> result_type;
	
	explicit divide_by(const Value2& v2) : v2(v2) {}
	
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
    typename tfunctor::divide_by<Value1, Value2, typename traits::algebraic_category<Value1>::type>::result_type
    inline divide_by(const Value1& value1, const Value2& value2)
    {
	return tfunctor::divide_by<Value1, Value2, typename traits::algebraic_category<Value1>::type>(value2)(value1);
    }
}

namespace vec {

    template <typename Value1, typename Value2>
    typename tfunctor::divide_by<Value1, Value2, typename traits::algebraic_category<Value1>::type>::result_type
    inline divide_by(const Value1& value1, const Value2& value2)
    {
	return tfunctor::divide_by<Value1, Value2, typename traits::algebraic_category<Value1>::type>(value2)(value1);
    }
}
	
} // namespace mtl

#endif // MTL_DIVIDE_BY_INCLUDE
