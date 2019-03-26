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

// With contributions from Cornelius Steinhardt

#ifndef MTL_MATRIX_EIGENVALUE_SYMMETRIC_INCLUDE
#define MTL_MATRIX_EIGENVALUE_SYMMETRIC_INCLUDE

#include <cmath>
#include <boost/utility.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/make_copy_or_reference.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/diagonal.hpp>
#include <boost/numeric/mtl/operation/givens.hpp>
#include <boost/numeric/mtl/operation/hessenberg.hpp>
#include <boost/numeric/mtl/operation/householder.hpp>
#include <boost/numeric/mtl/operation/qr.hpp>
#include <boost/numeric/mtl/operation/rank_one_update.hpp>
#include <boost/numeric/mtl/operation/signum.hpp>
#include <boost/numeric/mtl/operation/trans.hpp>

#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace mat {


/// Eigenvalues of symmetric matrix A with implicit QR algorithm
// Return Diagonalmatrix with eigenvalues as diag(A)
template <typename Matrix>
mtl::vec::dense_vector<typename Collection<Matrix>::value_type, mtl::vec::parameters<> >
inline qr_sym_imp(const Matrix& A)
{
    vampir_trace<5010> tracer;
    using std::abs; using mtl::signum; using mtl::real;
    typedef typename Collection<Matrix>::value_type   value_type;
    typedef typename Magnitude<value_type>::type      magnitude_type; // to multiply with 2 not 2+0i
    typedef typename Collection<Matrix>::size_type    size_type;
    size_type        ncols = num_cols(A), nrows = num_rows(A), N;
    value_type       zero= math::zero(A[0][0]), h00, h10, h11, beta, mu, a, b, tol;
    const magnitude_type two(2);
    Matrix           Q(nrows,ncols), H(nrows,ncols),  G(2,2);

    tol= 1.0e-8; // ????evtl ein Iterator wie bei den Gleichungssystemen, Problem: Keine Rechte Seite bzw b

    MTL_THROW_IF(ncols != nrows , matrix_not_square());

    // Hessenberg_form of Matrix A
    H= hessenberg(A);
    N= nrows;

    // QR_algo with implizit sym QR-step from Wilkinson
    while (1) {
	h00= H[N-2][N-2];
	h10= H[N-1][N-2];
	h11= H[N-1][N-1];

	//reduction, residuum and watch for breakdown
	if(abs(h10) < tol * abs(h11 + h00)) 
	    N--;	
	if (N < 2) 
	    break;
	
	// Wilkinson_shift
	beta= (h00 - h11) / two;   
	mu = h11 + (beta != zero ? beta - signum(beta) * sqrt(beta * beta + h10 * h10) : -h10);
	a= H[0][0] - mu, b= H[1][0];

	//implizit QR-step
	for (size_type k = 0; k < N - 1; k++) {
	    givens<Matrix>(H, a, b).trafo(k);
	    if (k < N - 2)
		a= H[k+1][k], b= H[k+2][k];	    
	}
    }
    return diagonal(H);
}


/// Evaluation of eigenvalues with QR-Algorithm of matrix A
// Return Diagonalmatrix with eigenvalues as diag(A)
template <typename Matrix>
mtl::vec::dense_vector<typename Collection<Matrix>::value_type, mtl::vec::parameters<> >
inline qr_algo(const Matrix& A, typename Collection<Matrix>::size_type itMax)
{
    vampir_trace<5011> tracer;
    typedef typename Collection<Matrix>::size_type    size_type;
    size_type        ncols = num_cols(A), nrows = num_rows(A);
    Matrix           Q(nrows, ncols), H(nrows, ncols), R(nrows, ncols);

    MTL_THROW_IF(ncols != nrows , matrix_not_square());

    H= hessenberg(A);
    for (size_type i = 0; i < itMax; i++) {
	boost::tie(Q, R)= qr_factors(H);
	H= R * Q;
    }
    return diagonal(H);
}


# ifdef MTL_SYMMETRIC_EIGENVALUE_WITH_QR

/// Calculation of eigenvalues of symmetric matrix A
template <typename Matrix>
mtl::vec::dense_vector<typename Collection<Matrix>::value_type, mtl::vec::parameters<> >
inline eigenvalue_symmetric(const Matrix& A, typename Collection<Matrix>::size_type itMax)
{
    return qr_algo(A, itMax == 0 ? num_rows(A) : itMax);
}

#else 

/// Calculation of eigenvalues of symmetric matrix A
template <typename Matrix>
mtl::vec::dense_vector<typename Collection<Matrix>::value_type, mtl::vec::parameters<> >
inline eigenvalue_symmetric(const Matrix& A, typename Collection<Matrix>::size_type)
{
    typedef dense2D<typename Collection<Matrix>::value_type, parameters<> >    arg_type;
    make_in_copy_or_reference<arg_type, Matrix>  copy_or_ref(A);
    return qr_sym_imp(copy_or_ref.value);
}

#endif

#if 0 // Too nasty to get it through all warnings :-!

/// Calculation of eigenvalues of symmetric matrix A
template <typename Matrix>
mtl::vec::dense_vector<typename Collection<Matrix>::value_type, mtl::vec::parameters<> >
inline eigenvalue_symmetric(const Matrix& A, 
			    typename Collection<Matrix>::size_type itMax= 0)
{
    vampir_trace<5012> tracer;
# ifdef MTL_SYMMETRIC_EIGENVALUE_WITH_QR
    return qr_algo(A, itMax == 0 ? num_rows(A) : itMax);
# else
    itMax= 0; // for not yelling at unused variable
    // qr_sym_imp works only with dense matrices of dynamic size, for other types copy
    typedef dense2D<typename Collection<Matrix>::value_type>    arg_type;
    make_in_copy_or_reference<arg_type, Matrix>  copy_or_ref(A);
    return qr_sym_imp(copy_or_ref.value);
# endif
}

#endif

template <typename Matrix>
mtl::vec::dense_vector<typename Collection<Matrix>::value_type, mtl::vec::parameters<> >
inline eigenvalue_symmetric(const Matrix& A)
{
    return eigenvalue_symmetric(A, 0);
} 


}} // namespace mtl::matrix


#endif // MTL_MATRIX_EIGENVALUE_SYMMETRIC_INCLUDE

