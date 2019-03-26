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

#ifndef MTL_MATRIX_MAKE_FAST_MULTI_VECTOR_EXPR_INCLUDE
#define MTL_MATRIX_MAKE_FAST_MULTI_VECTOR_EXPR_INCLUDE

#include <boost/numeric/mtl/utility/fast_multi_vector_expr.hpp>
#include <boost/numeric/mtl/matrix/mat_mat_plus_expr.hpp>
#include <boost/numeric/mtl/matrix/mat_mat_minus_expr.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/vec_vec_plus_expr.hpp>
#include <boost/numeric/mtl/vector/vec_vec_minus_expr.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl { namespace mat {


template <typename E1, typename E2>
typename mtl::traits::fast_multi_vector_expr< mv_mv_plus_expr<E1, E2> >::type
inline make_fast_multi_vector_expr(const mv_mv_plus_expr<E1, E2>& expr)
{
    typedef typename mtl::traits::fast_multi_vector_expr< mv_mv_plus_expr<E1, E2> >::type type;
    return type(make_fast_multi_vector_expr(expr.first), make_fast_multi_vector_expr(expr.second));
}

template <typename E1, typename E2>
typename mtl::traits::fast_multi_vector_expr< mv_mv_minus_expr<E1, E2> >::type
inline make_fast_multi_vector_expr(const mv_mv_minus_expr<E1, E2>& expr)
{
    typedef typename mtl::traits::fast_multi_vector_expr< mv_mv_minus_expr<E1, E2> >::type type;
    return type(make_fast_multi_vector_expr(expr.first), make_fast_multi_vector_expr(expr.second));
}

template <typename Functor, typename Matrix> 
typename mtl::traits::fast_multi_vector_expr< map_view<Functor, Matrix> >::type
inline make_fast_multi_vector_expr(const map_view<Functor, Matrix>& expr)
{
    typedef typename mtl::traits::fast_multi_vector_expr< map_view<Functor, Matrix> >::type type;
    return type(expr.functor, make_fast_multi_vector_expr(expr.ref));
}

template <typename Value1, typename Matrix>
typename mtl::traits::fast_multi_vector_expr< scaled_view<Value1, Matrix> >::type
inline make_fast_multi_vector_expr(const scaled_view<Value1, Matrix>& expr)
{
    typedef typename mtl::traits::fast_multi_vector_expr< scaled_view<Value1, Matrix> >::type type;
    return type(expr.functor.value(), make_fast_multi_vector_expr(expr.ref));
}

template <typename Value1, typename Matrix>
typename mtl::traits::fast_multi_vector_expr< rscaled_view<Value1, Matrix> >::type
inline make_fast_multi_vector_expr(const rscaled_view<Value1, Matrix>& expr)
{
    typedef typename mtl::traits::fast_multi_vector_expr< rscaled_view<Value1, Matrix> >::type type;
    return type(make_fast_multi_vector_expr(expr.ref), expr.functor.value());
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_MAKE_FAST_MULTI_VECTOR_EXPR_INCLUDE
