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

#ifndef MTL_MATRIX_LOWER_INCLUDE
#define MTL_MATRIX_LOWER_INCLUDE

namespace mtl { namespace mat {

namespace traits {

    template <typename Matrix>
    struct lower
    {
	typedef typename traits::bands<Matrix>::type type;
    };
}

/// Lower triangular matrix
template <typename Matrix> 
typename traits::lower<Matrix>::type
inline lower(const Matrix& A)
{
    return bands(A, std::numeric_limits<long>::min(), 1);
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_LOWER_INCLUDE
