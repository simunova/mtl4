// Copyright 2006. Peter Gottschling, Matthias Troyer, Rolf Bonderer

#ifndef LA_ETS_CONCEPTS_INCLUDE
#define LA_ETS_CONCEPTS_INCLUDE

#include <boost/numeric/linear_algebra/concepts.hpp>

#ifdef  __GXX_CONCEPTS__

//"Namespace ETS = Expression Templates Support":
// This file contains concepts that support the concepts defined in the namespace math.
// They are required since the math-concepts do not support ExpressionTemplates so far.
// The list of valid expressions is in fact infinite, so we just name some of them.
// Once ExpressionTemplates are supported by the math-concepts, the ets-concepts are no longer required.

namespace ets {

  auto concept Field<typename Element>
  {
    requires std::Assignable<Element, math::AdditivePartiallyInvertibleMonoid<Element>::unary_result_type>; //"x=-y" valid
  };
  
  auto concept VectorSpace<typename Vector, typename Scalar>
  {
    // valid expression: "vector2 += scalar*vector1"
    typename res_type_1;
    res_type_1 operator+=(Vector&, math::Multiplicable<Scalar, Vector>::result_type);
    
    // valid expression: "vector2 -= scalar*vector1"
    typename res_type_2;
    res_type_2 operator-=(Vector&, math::Multiplicable<Scalar, Vector>::result_type);
    
    //valid expression: "vector *= -scalar"
    typename res_type_3;
    res_type_3 operator*=(Vector&, const math::AdditivePartiallyInvertibleMonoid<Scalar>::unary_result_type&);
    
    //valid expression: "vector3 = vector1 + scalar*vector2"
    requires math::Addable<Vector, math::Multiplicable<Scalar, Vector>::result_type>; //"vector1+scalar*vector2" valid 
    requires std::Assignable<Vector, math::Addable<Vector, math::Multiplicable<Scalar, Vector>::result_type>::result_type>; //"vector3 = vector1 + scalar*vector2" valid
  
  };
   
} //namespace ets

#endif //__GXX_CONCEPTS__

#endif //LA_ETS_CONCEPTS_INCLUDE
