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

#ifndef MTL_TRAITS_IS_MULTI_VECTOR_EXPR_INCLUDE
#define MTL_TRAITS_IS_MULTI_VECTOR_EXPR_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_composable_vector.hpp>

namespace mtl { namespace traits {

/// Type trait whether an expression can be evaluated by multiple vector operations
template <typename T>
struct is_multi_vector_expr
  : boost::mpl::false_
{};

template <typename Vector>
struct is_multi_vector_expr< mtl::mat::multi_vector<Vector> >
  : boost::mpl::true_
{};

// template <typename E1, typename E2>
// struct is_multi_vector_expr< mtl::mat::mat_mat_asgn_expr<E1, E2> > 
//     : boost::mpl::bool_< is_multi_vector_expr<E1>::value && is_multi_vector_expr<E2>::value >
// {};

template <typename E1, typename E2>
struct is_multi_vector_expr< mtl::mat::mv_mv_plus_expr<E1, E2> > 
  : boost::mpl::bool_< is_multi_vector_expr<E1>::value && is_multi_vector_expr<E2>::value >
{};

template <typename E1, typename E2>
struct is_multi_vector_expr< mtl::mat::mat_mat_minus_expr<E1, E2> > 
    : boost::mpl::bool_< is_multi_vector_expr<E1>::value && is_multi_vector_expr<E2>::value >
{};

// template <typename E1, typename E2>
// struct is_multi_vector_expr< mtl::mat::mat_mat_ele_times_expr<E1, E2> > 
//     : boost::mpl::bool_< is_multi_vector_expr<E1>::value && is_multi_vector_expr<E2>::value >
// {};

template <typename Value1, typename Matrix>
struct is_multi_vector_expr< mtl::mat::scaled_view<Value1, Matrix> > 
  : is_multi_vector_expr<Matrix>
{};

template <typename Value1, typename Matrix>
struct is_multi_vector_expr< mtl::mat::rscaled_view<Value1, Matrix> > 
  : is_multi_vector_expr<Matrix>
{};

// ----------------------------------------------------


/// Type trait whether an expression can be evaluated by a single vector operation
template <typename T>
struct is_fast_multi_vector_expr
  : boost::mpl::false_
{};

template <typename Vector>
struct is_fast_multi_vector_expr< mtl::mat::multi_vector<Vector> >
  : is_composable_vector<Vector>
{};

// template <typename E1, typename E2>
// struct is_fast_multi_vector_expr< mtl::mat::mat_mat_asgn_expr<E1, E2> > 
//     : boost::mpl::bool_< is_fast_multi_vector_expr<E1>::value && is_fast_multi_vector_expr<E2>::value >
// {};

template <typename E1, typename E2>
struct is_fast_multi_vector_expr< mtl::mat::mv_mv_plus_expr<E1, E2> > 
  : boost::mpl::bool_< is_fast_multi_vector_expr<E1>::value && is_fast_multi_vector_expr<E2>::value >
{};

template <typename E1, typename E2>
struct is_fast_multi_vector_expr< mtl::mat::mat_mat_minus_expr<E1, E2> > 
    : boost::mpl::bool_< is_fast_multi_vector_expr<E1>::value && is_fast_multi_vector_expr<E2>::value >
{};

// template <typename E1, typename E2>
// struct is_fast_multi_vector_expr< mtl::mat::mat_mat_ele_times_expr<E1, E2> > 
//     : boost::mpl::bool_< is_fast_multi_vector_expr<E1>::value && is_fast_multi_vector_expr<E2>::value >
// {};

template <typename Value1, typename Matrix>
struct is_fast_multi_vector_expr< mtl::mat::scaled_view<Value1, Matrix> > 
  : is_fast_multi_vector_expr<Matrix>
{};

template <typename Value1, typename Matrix>
struct is_fast_multi_vector_expr< mtl::mat::rscaled_view<Value1, Matrix> > 
  : is_fast_multi_vector_expr<Matrix>
{};



}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_MULTI_VECTOR_EXPR_INCLUDE
