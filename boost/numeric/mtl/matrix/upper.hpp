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

#ifndef MTL_MATRIX_UPPER_INCLUDE
#define MTL_MATRIX_UPPER_INCLUDE

#include <boost/numeric/mtl/matrix/bands.hpp>
#include <limits>

namespace mtl { namespace mat {

namespace traits {

    template <typename Matrix>
    struct upper
    {
	typedef typename traits::bands<Matrix>::type type;
    };
}

/// Upper triangular matrix
template <typename Matrix> 
typename traits::upper<Matrix>::type
inline upper(const Matrix& A)
{
    return bands(A, 0, std::numeric_limits<long>::max());
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_UPPER_INCLUDE
