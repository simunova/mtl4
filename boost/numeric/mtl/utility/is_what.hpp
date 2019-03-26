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

#ifndef MTL_TRAITS_IS_WHAT_INCLUDE
#define MTL_TRAITS_IS_WHAT_INCLUDE

#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/bool.hpp>

#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/category.hpp>

namespace mtl { namespace traits {

// Matrix, vector and scalar deduced from ashape not from category
template <typename T>
struct is_scalar 
  : boost::is_same<typename ashape::ashape<T>::type, ashape::scal>
{};

template <typename Shape>
struct is_matrix_aux 
  : boost::mpl::false_
{};

template <typename Shape>
struct is_matrix_aux<ashape::mat<Shape> > 
  : boost::mpl::true_
{};

template <typename T>
struct is_matrix 
  : is_matrix_aux<typename ashape::ashape<T>::type>
{};

template <typename Shape>
struct is_vector_aux 
  : boost::mpl::false_
{};

template <typename Shape>
struct is_vector_aux<ashape::cvec<Shape> > 
  : boost::mpl::true_
{};

template <typename Shape>
struct is_vector_aux<ashape::rvec<Shape> > 
  : boost::mpl::true_
{};

template <typename T>
struct is_vector 
  : is_vector_aux<typename ashape::ashape<T>::type>
{};

template <typename T>
struct is_dense 
  : boost::is_base_of<tag::dense, typename category<T>::type> 
{};


template <typename T>
struct is_unevaluated 
  : boost::is_base_of<tag::unevaluated, typename category<T>::type> 
{};

template <typename T>
struct is_sparse 
  : boost::is_base_of<tag::sparse, typename category<T>::type> 
{};

// So far nothing is symmetric on the type level
template <typename T>
struct is_symmetric
    : boost::mpl::false_
{};


}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_WHAT_INCLUDE
