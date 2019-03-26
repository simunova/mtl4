// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_VECTOR_VEC_SCAL_MIXED_EXPR_INCLUDE
#define MTL_VECTOR_VEC_SCAL_MIXED_EXPR_INCLUDE

#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/operation/tfunctor_mixed.hpp>

namespace mtl { namespace vec {

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E1, map_view<tfunctor::left_plus<E1, typename E2::value_type>, E2> >::type
inline operator+(const E1& e1, const vec_expr<E2>& e2)
{
    typedef tfunctor::left_plus<E1, typename E2::value_type> ftype;
    typedef map_view<ftype, E2>                             type;

    return type(ftype(e1), static_cast<const E2&>(e2));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E2, map_view<tfunctor::right_plus<typename E1::value_type, E2>, E1> >::type
inline operator+(const vec_expr<E1>& e1, const E2& e2)
{
    typedef tfunctor::right_plus<typename E1::value_type, E2>ftype;
    typedef map_view<ftype, E1>                             type;

    return type(ftype(e2), static_cast<const E1&>(e1));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E1, map_view<tfunctor::left_minus<E1, typename E2::value_type>, E2> >::type
inline operator-(const E1& e1, const vec_expr<E2>& e2)
{
    typedef tfunctor::left_minus<E1, typename E2::value_type> ftype;
    typedef map_view<ftype, E2>                             type;

    return type(ftype(e1), static_cast<const E2&>(e2));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E2, map_view<tfunctor::right_minus<typename E1::value_type, E2>, E1> >::type
inline operator-(const vec_expr<E1>& e1, const E2& e2)
{
    typedef tfunctor::right_minus<typename E1::value_type, E2>ftype;
    typedef map_view<ftype, E1>                             type;

    return type(ftype(e2), static_cast<const E1&>(e1));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E1, map_view<tfunctor::left_min<E1, typename E2::value_type>, E2> >::type
inline min(const E1& e1, const vec_expr<E2>& e2)
{
    typedef tfunctor::left_min<E1, typename E2::value_type> ftype;
    typedef map_view<ftype, E2>                             type;

    return type(ftype(e1), static_cast<const E2&>(e2));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E2, map_view<tfunctor::right_min<typename E1::value_type, E2>, E1> >::type
inline min(const vec_expr<E1>& e1, const E2& e2)
{
    typedef tfunctor::right_min<typename E1::value_type, E2>ftype;
    typedef map_view<ftype, E1>                             type;

    return type(ftype(e2), static_cast<const E1&>(e1));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E1, map_view<tfunctor::left_max<E1, typename E2::value_type>, E2> >::type
inline max(const E1& e1, const vec_expr<E2>& e2)
{
    typedef tfunctor::left_max<E1, typename E2::value_type> ftype;
    typedef map_view<ftype, E2>                             type;

    return type(ftype(e1), static_cast<const E2&>(e2));
}

template <typename E1, typename E2>
typename mtl::traits::enable_if_scalar<E2, map_view<tfunctor::right_max<typename E1::value_type, E2>, E1> >::type
inline max(const vec_expr<E1>& e1, const E2& e2)
{
    typedef tfunctor::right_max<typename E1::value_type, E2>ftype;
    typedef map_view<ftype, E1>                             type;

    return type(ftype(e2), static_cast<const E1&>(e1));
}

}} // namespace mtl::vector

#endif // MTL_VECTOR_VEC_SCAL_MIXED_EXPR_INCLUDE
