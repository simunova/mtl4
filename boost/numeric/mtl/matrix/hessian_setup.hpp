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

#ifndef MTL_HESSIAN_MATRIX_UTILITIES_INCLUDE
#define MTL_HESSIAN_MATRIX_UTILITIES_INCLUDE

#include <cmath>
#include <iostream>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/operation/entry_similar.hpp>

namespace mtl { namespace mat {

/// Fills a matrix A with A[i][j] = factor * (i + j)
/** Intended for dense matrices.
    Works on sparse matrices with inserter but is very expensive. **/
template <typename Matrix, typename Value>
void hessian_setup(Matrix& A, Value factor)
{
    typedef typename Collection<Matrix>::value_type    value_type;
    typedef typename Collection<Matrix>::size_type     size_type;
    { // extra block unfortunately needed for VS2013
	inserter<Matrix> ins(A, num_cols(A));

	for (size_type r = 0; r < num_rows(A); r++)
	for (size_type c = 0; c < num_cols(A); c++)
	    ins[r][c] << factor * (value_type(r) + value_type(c));
    }
}

namespace impl {

    /*
    - Check matrix product C = A * B with:
      - A is MxN, B is NxL, C is MxL
      - with matrices a_ij = i+j, b_ij = 2(i+j); 
      - c_ij = 1/3 N (1 - 3i - 3j + 6ij - 3N + 3iN + 3jN + 2N^2).

    */
    // Not really generic
    template <typename Value>
    double inline hessian_product_i_j (Value i, Value j, Value N)
    {
        return 1.0/3.0 * N * (1.0 - 3*i - 3*j + 6*i*j - 3*N + 3*i*N + 3*j*N + 2*N*N);
    }
        
    template <typename Value>
    inline bool similar_values(Value x, Value y) 
    {
        using std::abs; using std::max;
        return abs(x - y) / max(abs(x), abs(y)) < 0.000001;
    }

    template <typename Matrix>
    void inline check_entry(Matrix const& C, std::size_t r, std::size_t c, 
			    std::size_t reduced_dim, double factor)
    {
	if (!entry_similar(C, r, c, factor * hessian_product_i_j(r, c, reduced_dim), 0.00001)) {
	    std::cerr << "Result in C[" << r << "][" << c << "] should be " 
		      << factor * hessian_product_i_j(r, c, reduced_dim)
		      << " but is " << C[r][c] << "\n";
	    MTL_THROW(unexpected_result());
	}
    }

} // impl       


/// Check if matrix C is A * B with A and B set by hessian_setup
/** C has dimensions M x L and reduced_dim is N, see hessian_setup. **/
template <typename Matrix>
void check_hessian_matrix_product(Matrix const& C, std::size_t reduced_dim, double factor= 1.0)
{
    if (num_rows(C) * num_cols(C) == 0) return; // otherwise out of range

    impl::check_entry(C, 0, 0, reduced_dim, factor);
    impl::check_entry(C, 0, num_cols(C)-1, reduced_dim, factor);
    impl::check_entry(C, num_rows(C)-1, 0, reduced_dim, factor);
    impl::check_entry(C, num_rows(C)-1, num_cols(C)-1, reduced_dim, factor);
    impl::check_entry(C, num_rows(C)/2, num_cols(C)/2, reduced_dim, factor);
}

} // namespace matrix;

using mat::hessian_setup;
using mat::check_hessian_matrix_product;

} // namespace mtl

#endif // MTL_HESSIAN_MATRIX_UTILITIES_INCLUDE
