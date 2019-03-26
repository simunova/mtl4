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

#ifndef MTL_MATRIX_INV_INCLUDE
#define MTL_MATRIX_INV_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/identity.hpp>
#include <boost/numeric/mtl/operation/upper_trisolve.hpp>
#include <boost/numeric/mtl/operation/lu.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/unit_vector.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

namespace traits {

    /// Return type of inv(Matrix)
    /** Might be specialized later for the sake of efficiency **/
    template <typename Matrix>
    struct inv
    {
	typedef typename Collection<Matrix>::value_type    value_type;
	typedef ::mtl::mat::dense2D<value_type>         type;
    };
	
} // traits

/// Invert upper triangular matrix
template <typename Matrix, typename MatrixOut>
void inv_upper(Matrix const& A, MatrixOut& Inv)
{
    vampir_trace<5019> tracer;
    typedef typename Collection<Matrix>::value_type    value_type;
    typedef typename Collection<Matrix>::size_type     size_type;
   
    const size_type N= num_rows(A);
    MTL_DEBUG_THROW_IF(num_cols(A) != N, matrix_not_square());
    MTL_DEBUG_THROW_IF(N != num_rows(Inv) || num_cols(A) != num_cols(Inv), incompatible_size());

    Inv= math::zero(value_type());

    for (size_type k= 0; k < N; ++k) {
	irange r(k+1);
	typename mtl::ColumnInMatrix<MatrixOut>::type col_k(Inv[r][k]);
	upper_trisolve(A[r][r], vec::unit_vector<value_type>(k, k+1), col_k, mtl::tag::regular_diagonal());
    }
}


/// Invert upper triangular matrix
template <typename Matrix>
inline typename traits::inv<Matrix>::type 
inv_upper(Matrix const& A)
{
    typedef typename Collection<Matrix>::size_type     size_type;
    const size_type N= num_rows(A);
    typename traits::inv<Matrix>::type Inv(N, N);
    inv_upper(A, Inv);
    return Inv;
}

#if 0
/// Invert lower triangular matrix
template <typename Matrix, typename MatrixOut>
inline void inv_lower(Matrix const& A, MatrixOut& Inv)
{
    vampir_trace<5020> tracer;
    typedef typename Collection<Matrix>::value_type    value_type;
    typedef typename Collection<Matrix>::size_type     size_type;

    const size_type N= num_rows(A);
    MTL_DEBUG_THROW_IF(num_cols(A) != N, matrix_not_square());
    MTL_DEBUG_THROW_IF(N != num_rows(Inv) || num_cols(A) != num_cols(Inv), incompatible_size());

    Inv= math::zero(value_type());

    for (size_type k= 0; k < N; ++k) {
	irange r(k, N);
	typename mtl::ColumnInMatrix<MatrixOut>::type col_k(Inv[r][k]);
	lower_trisolve(A[r][r], unit_vector<value_type>(0, N-k), col_k, mtl::tag::regular_diagonal());
    }
}

template <typename Matrix>
typename traits::inv<Matrix>::type
inline inv_lower(Matrix const& A)
{
    typedef typename Collection<Matrix>::size_type     size_type;
    const size_type N= num_rows(A);
    typename traits::inv<Matrix>::type Inv(N, N);
    inv_lower(A, Inv);
    return Inv;
}
#endif

#if 1
/// Invert lower triangular matrix
template <typename Matrix>
typename traits::inv<Matrix>::type
inline inv_lower(Matrix const& A)
{
    vampir_trace<5020> tracer;
    Matrix T(trans(A)); // Shouldn't be needed
    return typename traits::inv<Matrix>::type(trans(inv_upper(T)));
}
#endif




/// Invert matrix
/** Uses pivoting LU factorization and triangular inversion
    \sa \ref lu, \ref inv_upper, \ref inv_lower **/
template <typename Matrix, typename MatrixOut>
inline void inv(Matrix const& A, MatrixOut& Inv)
{
    vampir_trace<5021> tracer;
    typedef typename Collection<Matrix>::size_type     size_type;
    typedef typename Collection<Matrix>::value_type    value_type;
    typedef typename traits::inv<Matrix>::type         result_type;

    const size_type N= num_rows(A);
    MTL_THROW_IF(num_cols(A) != num_cols(A), matrix_not_square());
    MTL_DEBUG_THROW_IF(N != num_rows(Inv) || num_cols(A) != num_cols(Inv), incompatible_size());

    if (N == 1) {
	Inv[0][0]= value_type(1) / A[0][0];
	return;
    }

    result_type                    PLU(A);
    mtl::dense_vector<size_type, vec::parameters<> >   Pv(num_rows(A));

    lu(PLU, Pv);
    result_type  PU(upper(PLU)), PL(strict_lower(PLU));
    for (size_type i= 0; i < num_rows(A); i++)
	PL[i][i]= value_type(1);

    Inv= inv_upper(PU) * inv_lower(PL) * permutation(Pv);
}

template <typename Matrix>
typename traits::inv<Matrix>::type
inline inv(Matrix const& A)
{
    typedef typename Collection<Matrix>::size_type     size_type;
    const size_type N= num_rows(A);
    typename traits::inv<Matrix>::type Inv(N, N);
    inv(A, Inv);
    return Inv;
}

}} // namespace mtl::matrix

#endif // MTL_MATRIX_INV_INCLUDE
