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

#ifndef MTL_TRAITS_TRANSPOSED_ORIENTATION_INCLUDE
#define MTL_TRAITS_TRANSPOSED_ORIENTATION_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Orientation type for transposed matrices and vectors
template <class T> struct transposed_orientation {};

template<> struct transposed_orientation<tag::row_major> 
{
    typedef tag::col_major type; 
};

template<> struct transposed_orientation<tag::col_major> 
{
    typedef tag::row_major type; 
};


}} // namespace mtl::traits

#endif // MTL_TRAITS_TRANSPOSED_ORIENTATION_INCLUDE
