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

#ifndef MTL_TRAITS_EVAL_DENSE_INCLUDE
#define MTL_TRAITS_EVAL_DENSE_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {


template <typename T>
struct eval_dense
    : boost::mpl::false_
{};

template <typename Value, typename Parameter>
struct eval_dense< mtl::vec::dense_vector<Value, Parameter> >
    : boost::mpl::true_
{};

template <typename Value, typename Parameter>
struct eval_dense< mtl::mat::dense2D<Value, Parameter> >
    : boost::mpl::true_
{};

template <typename Value, std::size_t Mask, typename Parameter>
struct eval_dense< mtl::mat::morton_dense<Value, Mask, Parameter> >
    : boost::mpl::true_
{};



template <typename Value1, typename Vector>
struct eval_dense< mtl::vec::scaled_view<Value1, Vector> > 
    : eval_dense<Vector>
{};

template <typename Value1, typename Vector>
struct eval_dense< mtl::vec::rscaled_view<Value1, Vector> > 
    : eval_dense<Vector>
{};



template <typename E1, typename E2>
struct eval_dense< mtl::mat::mat_mat_asgn_expr<E1, E2> > 
    : boost::mpl::bool_< eval_dense<E1>::value && eval_dense<E2>::value >
{};

template <typename E1, typename E2>
struct eval_dense< mtl::mat::mat_mat_plus_expr<E1, E2> > 
    : boost::mpl::bool_< eval_dense<E1>::value && eval_dense<E2>::value >
{};

template <typename E1, typename E2>
struct eval_dense< mtl::mat::mat_mat_minus_expr<E1, E2> > 
    : boost::mpl::bool_< eval_dense<E1>::value && eval_dense<E2>::value >
{};

template <typename E1, typename E2>
struct eval_dense< mtl::mat::mat_mat_ele_times_expr<E1, E2> > 
    : boost::mpl::bool_< eval_dense<E1>::value && eval_dense<E2>::value >
{};

template <typename Value1, typename Matrix>
struct eval_dense< mtl::mat::scaled_view<Value1, Matrix> > 
    : eval_dense<Matrix>
{};

template <typename Matrix>
struct eval_dense< mtl::mat::negate_view<Matrix > > 
    : eval_dense<Matrix>
{};

template <typename Matrix>
struct eval_dense< mtl::mat::real_view<Matrix > > 
    : eval_dense<Matrix>
{};

template <typename Matrix>
struct eval_dense< mtl::mat::imag_view<Matrix > > 
    : eval_dense<Matrix>
{};

template <typename Functor, typename Matrix>
struct eval_dense< mtl::mat::map_view<Functor, Matrix > > 
    : eval_dense<Matrix>
{};

template <typename Value1, typename Matrix>
struct eval_dense< mtl::mat::rscaled_view<Value1, Matrix> > 
    : eval_dense<Matrix>
{};


}} // namespace mtl::traits

#endif // MTL_TRAITS_EVAL_DENSE_INCLUDE
