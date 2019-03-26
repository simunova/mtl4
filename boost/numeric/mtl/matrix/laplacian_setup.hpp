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

#ifndef MTL_LAPLACIAN_SETUP_INCLUDE
#define MTL_LAPLACIAN_SETUP_INCLUDE

#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

/// Setup a matrix according to a Laplacian equation on a 2D-grid using a five-point-stencil
/** Intended for sparse matrices but works also with dense matrices. Changes the size of
    the matrix \f$m\cdot n\times m\cdot n\f$. **/
template <typename Matrix>
inline void laplacian_setup(Matrix& A, unsigned m, unsigned n)
{
    vampir_trace<3063> tracer;
    A.change_dim(m*n, m*n);
    set_to_zero(A);
    { // extra block unfortunately needed for VS2013
	inserter<Matrix>      ins(A, 5);
	for (unsigned i = 0; i < m; i++)
	for (unsigned j = 0; j < n; j++) {
	    typename Collection<Matrix>::value_type four(4.0), minus_one(-1.0);
	    unsigned row = i * n + j;
	    ins(row, row) << four;
	    if (j < n - 1) ins(row, row + 1) << minus_one;
	    if (i < m - 1) ins(row, row + n) << minus_one;
	    if (j > 0) ins(row, row - 1) << minus_one;
	    if (i > 0) ins(row, row - n) << minus_one;
	}
    }
}

}} // namespace mtl::matrix

#endif // MTL_LAPLACIAN_SETUP_INCLUDE
