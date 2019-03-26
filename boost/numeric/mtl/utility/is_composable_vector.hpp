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

#ifndef MTL_TRAITS_IS_COMPOSABLE_VECTOR_INCLUDE
#define MTL_TRAITS_IS_COMPOSABLE_VECTOR_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

template <typename T>
struct is_composable_vector
  : boost::mpl::false_
{};

template <typename Value, typename Parameters>
struct is_composable_vector< mtl::vec::dense_vector<Value, Parameters> >
  : boost::mpl::true_
{};

	

}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_COMPOSABLE_VECTOR_INCLUDE
