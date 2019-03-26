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

#ifndef MTL_TRAITS_POS_TYPE_INCLUDE
#define MTL_TRAITS_POS_TYPE_INCLUDE

#include <utility>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>

namespace mtl { namespace traits {

// Not a matrix -> should be a vector
template <typename T, bool IsMatrix>
struct pos_type_aux
{
    typedef typename Collection<T>::size_type type;
};

// Matrix specialization
template <typename T>
struct pos_type_aux<T, true>
{
  private:
    typedef typename Collection<T>::size_type size_type;
  public:
    typedef std::pair<size_type, size_type>   type;
};

/// Type trait for position type
/** This is size_type for a vector and for a matrix a pair of size_types.
    N.B.: The implementation considers everything that is not a matrix as vector, to be improved one day. **/
template <typename T>
struct pos_type
  : pos_type_aux<T, is_matrix<T>::value>
{};


}} // namespace mtl::traits

#endif // MTL_TRAITS_POS_TYPE_INCLUDE
