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

#ifndef MTL_TRAITS_FAST_MULTI_VECTOR_EXPR_INCLUDE
#define MTL_TRAITS_FAST_MULTI_VECTOR_EXPR_INCLUDE

namespace mtl { namespace traits {


/// Type trait whether an expression can be evaluated by a single vector operation
template <typename T>
struct fast_multi_vector_expr
{};

template <typename Vector>
struct fast_multi_vector_expr< mtl::mat::multi_vector<Vector> >
{
    typedef Vector type;
};

// template <typename Vector>
// struct fast_multi_vector_expr< const mtl::mat::multi_vector<Vector> >
// {
//     typedef const Vector type;
// };

// template <typename E1, typename E2>
// struct fast_multi_vector_expr< mtl::mat::mat_mat_asgn_expr<E1, E2> > 
//     : boost::mpl::bool_< fast_multi_vector_expr<E1>::value && fast_multi_vector_expr<E2>::value >
// {};

template <typename E1, typename E2>
struct fast_multi_vector_expr< mtl::mat::mv_mv_plus_expr<E1, E2> > 
{
    typedef typename fast_multi_vector_expr<E1>::type V1;
    typedef typename fast_multi_vector_expr<E2>::type V2;
    typedef typename mtl::sfunctor::plus<typename V1::value_type, typename V2::value_type> Functor;

    typedef mtl::vec::vec_vec_pmop_expr<V1, V2, Functor> type;
};

template <typename E1, typename E2>
struct fast_multi_vector_expr< mtl::mat::mat_mat_minus_expr<E1, E2> > 
{
    typedef typename fast_multi_vector_expr<E1>::type V1;
    typedef typename fast_multi_vector_expr<E2>::type V2;
    typedef typename mtl::sfunctor::minus<typename V1::value_type, typename V2::value_type> Functor;

    typedef mtl::vec::vec_vec_pmop_expr<V1, V2, Functor> type;
};

// template <typename E1, typename E2>
// struct fast_multi_vector_expr< mtl::mat::mat_mat_ele_times_expr<E1, E2> > 
// {
//     typedef mtl::mat::mv_mv_ele_times_expr<typename fast_multi_vector_expr<E1>::type,
// 					      typename fast_multi_vector_expr<E2>::type> type;
// };

template <typename Functor, typename Matrix> 
struct fast_multi_vector_expr< mtl::mat::map_view<Functor, Matrix> >
{
    typedef mtl::vec::map_view<Functor, typename fast_multi_vector_expr<Matrix>::type> type;
};


template <typename Value1, typename Matrix>
struct fast_multi_vector_expr< mtl::mat::scaled_view<Value1, Matrix> > 
{
    typedef mtl::vec::scaled_view<Value1, typename fast_multi_vector_expr<Matrix>::type> type;
};

template <typename Value1, typename Matrix>
struct fast_multi_vector_expr< mtl::mat::rscaled_view<Value1, Matrix> > 
{
    typedef mtl::vec::rscaled_view<Value1, typename fast_multi_vector_expr<Matrix>::type> type;
};


}} // namespace mtl::traits

#endif // MTL_TRAITS_FAST_MULTI_VECTOR_EXPR_INCLUDE
