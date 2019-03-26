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

#ifndef MTL_TRAITS_ROOT_INCLUDE
#define MTL_TRAITS_ROOT_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>

namespace mtl { namespace traits {


/// Type trait to reduce types to their essentials by removing const, reference, ... and gearing derived types to their bases
template <typename T>
struct root
{
    typedef T        type;
};

// ==========================
// Remove language attributes
// ==========================

template <typename T>
struct root<T&>
  : public root<T> {};

template <typename T>
struct root<const T>
  : public root<T> {};

// Redundant specialization to make xlc++ happy
template <typename T, int R, int C>
struct root<const T[R][C]>
  : public root<T[R][C]> {};


// ============
// Base classes
// ============

// Implicit dense matrices

template <typename Value>
struct root<mtl::mat::ones_matrix<Value> >
{
    typedef mtl::mat::implicit_dense<mtl::mat::ones_functor<Value> > type;
};

template <typename Value>
struct root<mtl::mat::hilbert_matrix<Value> >
{
    typedef mtl::mat::implicit_dense<mtl::mat::hilbert_functor<Value> > type;
};

template <typename Vector1, typename Vector2>
struct root<mtl::mat::outer_product_matrix<Vector1, Vector2> >
{
    typedef mtl::mat::implicit_dense<mtl::mat::outer_product_functor<Vector1, Vector2> > type;
};

// Matrix map views

template <typename Scaling, typename Matrix>
struct root< mtl::mat::scaled_view<Scaling, Matrix> >
{
    typedef mtl::mat::map_view<tfunctor::scale<Scaling, typename Matrix::value_type>, Matrix> type;
};

template <typename Matrix>
struct root< mtl::mat::conj_view<Matrix> >
{
    typedef mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, Matrix> type;
};

template <typename Matrix>
struct root< mtl::mat::negate_view<Matrix> >
{
    typedef mtl::mat::map_view<sfunctor::negate<typename Matrix::value_type>, Matrix> type;
};

template <typename Matrix>
struct root< mtl::mat::imag_view<Matrix> >
{
    typedef mtl::mat::map_view<sfunctor::imag<typename Matrix::value_type>, Matrix> type;
};

template <typename Matrix>
struct root< mtl::mat::real_view<Matrix> >
{
    typedef mtl::mat::map_view<sfunctor::real<typename Matrix::value_type>, Matrix> type;
};

template <typename Matrix>
struct root< mtl::mat::exp_view<Matrix> >
{
    typedef mtl::mat::map_view<sfunctor::exp<typename Matrix::value_type>, Matrix> type;
};

template <typename Matrix>
struct root< mtl::mat::hermitian_view<Matrix> >
{
    typedef mtl::mat::map_view<sfunctor::conj<typename Matrix::value_type>, mtl::mat::transposed_view<Matrix> > type;
};

template <typename Matrix, typename RScaling>
struct root< mtl::mat::rscaled_view<Matrix, RScaling> >
{
    typedef mtl::mat::map_view<tfunctor::rscale<typename Matrix::value_type, RScaling>, Matrix> type;
};

template <typename Matrix, typename Divisor>
struct root< mtl::mat::divide_by_view<Matrix, Divisor> >
{
    typedef mtl::mat::map_view<tfunctor::divide_by<typename Matrix::value_type, Divisor>, Matrix> type;
};

// Matrix operations
template <typename M1, typename M2> 
struct root< mtl::mat::mat_mat_plus_expr<M1, M2> >
{
    typedef mtl::sfunctor::plus<typename Collection<M1>::value_type, typename Collection<M2>::value_type> f_type;
    typedef mtl::mat::mat_mat_op_expr<M1, M2, f_type> type;
};

template <typename M1, typename M2> 
struct root< mtl::mat::mv_mv_plus_expr<M1, M2> >
{
    typedef typename root< mtl::mat::mat_mat_plus_expr<M1, M2> >::type type;
};

template <typename M1, typename M2> 
struct root< mtl::mat::mat_mat_minus_expr<M1, M2> >
{
    typedef mtl::sfunctor::minus<typename Collection<M1>::value_type, typename Collection<M2>::value_type> f_type;
    typedef mtl::mat::mat_mat_op_expr<M1, M2, f_type> type;
};

template <typename M1, typename M2> 
struct root< mtl::mat::mv_mv_minus_expr<M1, M2> >
{
    typedef typename root< mtl::mat::mat_mat_minus_expr<M1, M2> >::type type;
};

template <typename M1, typename M2> 
struct root< mtl::mat::mat_mat_times_expr<M1, M2> >
{
    typedef mtl::sfunctor::times<typename Collection<M1>::value_type, typename Collection<M2>::value_type> f_type;
    typedef mtl::mat::mat_mat_op_expr<M1, M2, f_type> type;
};

template <typename M1, typename M2> 
struct root< mtl::mat::mat_mat_ele_times_expr<M1, M2> >
{
    typedef mtl::sfunctor::times<typename Collection<M1>::value_type, typename Collection<M2>::value_type> f_type;
    typedef mtl::mat::mat_mat_op_expr<M1, M2, f_type> type;
};


// Vector assignment expressions

template <typename E1, typename E2>
struct root< vec::vec_vec_asgn_expr<E1, E2> >
{
    typedef vec::vec_vec_aop_expr< E1, E2, mtl::sfunctor::assign<typename E1::value_type, typename E2::value_type> > type;
};

template <typename E1, typename E2>
struct root< vec::vec_vec_plus_asgn_expr<E1, E2> >
{
    typedef vec::vec_vec_aop_expr< E1, E2, mtl::sfunctor::plus_assign<typename E1::value_type, typename E2::value_type> > type;
};

template <typename E1, typename E2>
struct root< vec::vec_vec_minus_asgn_expr<E1, E2> >
{
    typedef vec::vec_vec_aop_expr< E1, E2, mtl::sfunctor::minus_assign<typename E1::value_type, typename E2::value_type> > type;
};

template <typename E1, typename E2>
struct root< vec::vec_scal_asgn_expr<E1, E2> >
{
    typedef vec::vec_scal_aop_expr< E1, E2, mtl::sfunctor::assign<typename E1::value_type, E2> > type;
};

template <typename E1, typename E2>
struct root< vec::vec_scal_times_asgn_expr<E1, E2> >
{
    typedef vec::vec_scal_aop_expr< E1, E2, mtl::sfunctor::times_assign<typename E1::value_type, E2> > type;
};

template <typename E1, typename E2>
struct root< vec::vec_scal_div_asgn_expr<E1, E2> >
{
    typedef vec::vec_scal_aop_expr< E1, E2, mtl::sfunctor::divide_assign<typename E1::value_type, E2> > type;
};

template <typename Scaling, typename Vector>
struct root< vec::scaled_view<Scaling, Vector> >
{
    typedef vec::map_view<tfunctor::scale<Scaling, typename Vector::value_type>, Vector> type;
};

template <typename Vector, typename RScaling>
struct root< vec::rscaled_view<Vector, RScaling> >
{
    typedef vec::map_view<tfunctor::rscale<typename Vector::value_type, RScaling>, Vector> type;
};

template <typename Vector, typename Divisor>
struct root< vec::divide_by_view<Vector, Divisor> >
{
    typedef vec::map_view<tfunctor::divide_by<typename Vector::value_type, Divisor>, Vector> type;
};

template <typename Vector>
struct root< vec::conj_view<Vector> >
{
    typedef vec::map_view<mtl::sfunctor::conj<typename Vector::value_type>, Vector> type;
};

template <typename Vector>
struct root< vec::negate_view<Vector> >
{
    typedef vec::map_view<mtl::sfunctor::negate<typename Vector::value_type>, Vector> type;
};

template <unsigned BSize, typename Vector>
struct root< vec::unrolled1<BSize, Vector> >
{
    typedef Vector type;
};



#if 0 // template
struct root
{
    typedef  type;
};
#endif


}} // namespace mtl::traits

#endif // MTL_TRAITS_ROOT_INCLUDE
