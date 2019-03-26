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

#ifndef MTL_TRAITS_VIEWED_COLLECTION_INCLUDE
#define MTL_TRAITS_VIEWED_COLLECTION_INCLUDE

#include <boost/type_traits/remove_const.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>

namespace mtl { namespace traits {

/// Type of the viewed matrix or vector, by default itself
template <typename T>
struct viewed_collection
{
    MTL_STATIC_ASSERT(is_matrix<T>::value, "Currently only matrices are supported.");
    typedef typename boost::remove_const<T>::type type;
};

template <typename Matrix>
struct viewed_collection< mat::conj_view<Matrix> >
  : viewed_collection<Matrix>
{};

template <typename Matrix>
struct viewed_collection< mat::transposed_view<Matrix> >
  : viewed_collection<Matrix>
{};

template <typename Matrix>
struct viewed_collection< mat::hermitian_view<Matrix> >
  : viewed_collection<Matrix>
{};

// add vector stuff

}} // namespace mtl::traits

#endif // MTL_TRAITS_VIEWED_COLLECTION_INCLUDE
