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

#ifndef MTL_MATRIX_STRICT_UPPER_INCLUDE
#define MTL_MATRIX_STRICT_UPPER_INCLUDE

#include <boost/numeric/mtl/matrix/bands.hpp>

namespace mtl { namespace mat {

namespace traits {

    template <typename Matrix>
    struct strict_upper
    {
	typedef typename bands<Matrix>::type type;
    };
}

///  Strict upper triangle matrix
template <typename Matrix> 
typename traits::strict_upper<Matrix>::type
inline strict_upper(const Matrix& A)
{
    return bands(A, 1, std::numeric_limits<long>::max());
}

/// Triangle-upper starting at off-diagonoal \p d (for compatibility with matlib)
template <typename Matrix> 
typename traits::strict_upper<Matrix>::type
inline triu(const Matrix& A, long d= 0)
{
    return bands(A, d, std::numeric_limits<long>::max());
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_STRICT_UPPER_INCLUDE
