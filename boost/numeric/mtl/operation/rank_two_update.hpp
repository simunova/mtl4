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

#ifndef MTL_RANK_TWO_UPDATE_INCLUDE
#define MTL_RANK_TWO_UPDATE_INCLUDE

#include <boost/numeric/mtl/operation/rank_one_update.hpp>

namespace mtl { namespace mat {

/// Rank-two update: rank_two_update(A, x, y) computes A+= x * conj(y)^T + y * conj(x)^T
/** The current implementation works for column and row vectors (although
    the notation above refers to column vectors). **/
template <typename Matrix, typename VectorX, typename VectorY>
inline void rank_two_update(Matrix& matrix, const VectorX& x, const VectorY& y)
{
    rank_one_update(matrix, x, y);
    rank_one_update(matrix, y, x);
}

}} // namespace mtl::matrix

#endif // MTL_RANK_TWO_UPDATE_INCLUDE
