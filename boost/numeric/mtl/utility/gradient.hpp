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

#ifndef MTL_TRAITS_GRADIENT_INCLUDE
#define MTL_TRAITS_GRADIENT_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/linear_operator.hpp>

namespace mtl { namespace traits {

/// Type trait returning type for gradient of function with T as argument
template <typename T>
struct gradient 
  : boost::mpl::if_<is_scalar<T>, mtl::dense_vector<T>, tag::unknown>
{};	

template <typename Value, typename Para>
struct gradient<mtl::dense_vector<Value, Para> >
  : linear_operator<mtl::dense_vector<Value, Para>, mtl::dense_vector<Value, Para> >
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_GRADIENT_INCLUDE
