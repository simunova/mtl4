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

#ifndef MTL_MATRIX_IDENTITY_INCLUDE
#define MTL_MATRIX_IDENTITY_INCLUDE

// #include <boost/numeric/linear_algebra/identity.hpp>
// #include <boost/numeric/mtl/mtl_fwd.hpp>
// #include <boost/numeric/mtl/matrix/parameter.hpp>
// #include <boost/numeric/mtl/matrix/diagonal_setup.hpp>

#include <boost/numeric/mtl/matrix/identity2D.hpp>

namespace mtl { namespace mat {

inline identity2D identity(std::size_t nrows, std::size_t ncols)
{
    return identity2D(nrows, ncols);
}


inline identity2D identity(std::size_t nrows)
{
    return identity2D(nrows, nrows);
}

}} // namespace mtl::matrix

#endif // MTL_MATRIX_IDENTITY_INCLUDE
