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

#ifndef MTL_TRAITS_IS_VECTOR_REDUCTION_INCLUDE
#define MTL_TRAITS_IS_VECTOR_REDUCTION_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/mpl/bool.hpp>

namespace mtl { namespace traits {

/// Type trait to check whether \p T is a vector reduction (i.e. an ET class for its lazy evaluation)
template <typename T>
struct is_vector_reduction : boost::mpl::false_ {};

template <unsigned long Unroll, typename Vector1, typename Vector2, typename ConjOpt>
struct is_vector_reduction<mtl::vec::dot_class<Unroll, Vector1, Vector2, ConjOpt> >
  : boost::mpl::true_ {};

template<typename Vector, typename Functor>
struct is_vector_reduction<mtl::vec::lazy_reduction<Vector, Functor> >
  : boost::mpl::true_ {};



}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_VECTOR_REDUCTION_INCLUDE
