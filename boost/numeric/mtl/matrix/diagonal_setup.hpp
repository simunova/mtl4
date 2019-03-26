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

#ifndef MTL_DIAGONAL_SETUP_INCLUDE
#define MTL_DIAGONAL_SETUP_INCLUDE

#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

namespace mtl { namespace mat {

/// Setup a matrix to a multiple of the unity matrix
/** Intended for sparse matrices but works also with dense matrices. 
    If the value is 0 the matrix is only zeroed out, whereby
    a sparse matrix will be empty after this operation,
    i.e. the zeros on the diagonal are not explicitly stored.
    The diagonal in its generalized form is the set of entries with equal row and column
    index (since r6843, older revision considered it erroneous to store
    a non-zero scalar to a non-square matrix).
 **/
template <typename Matrix, typename Value>
inline void diagonal_setup(Matrix& matrix, const Value& value)
{
    using std::min;
    if (num_rows(matrix) == 0 || num_cols(matrix) == 0) 
	return;

    set_to_zero(matrix);
    inserter<Matrix>      ins(matrix, 1);
    for (typename Collection<Matrix>::size_type i= 0, n= min(num_rows(matrix), num_cols(matrix)); i < n; ++i)
	ins[i][i] << value;
}

}} // namespace mtl::matrix

#endif // MTL_DIAGONAL_SETUP_INCLUDE
