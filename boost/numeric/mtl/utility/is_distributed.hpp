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

#ifndef MTL_TRAITS_IS_DISTRIBUTED_INCLUDE
#define MTL_TRAITS_IS_DISTRIBUTED_INCLUDE

namespace mtl { namespace traits {

template <typename T> 
struct is_distributed_aux 
  : boost::mpl::false_       // by default false
{};

/// Meta-function whether a certain type is distributed
template <typename T> 
struct is_distributed 
  : is_distributed_aux<typename root<T>::type>
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_DISTRIBUTED_INCLUDE
