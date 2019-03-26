// Software License for MTL
//
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// With contributions from Cornelius Steinhardt
//
// This file is part of the Matrix Template Library
//
// See also license.mtl.txt in the distribution.


#ifndef MTL_MATRIX_QR_INCLUDE
#define MTL_MATRIX_QR_INCLUDE

#include <cmath>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/operation/householder.hpp>
#include <boost/numeric/mtl/operation/rank_one_update.hpp>
#include <boost/numeric/mtl/operation/trans.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {


/// QR-Factorization of matrix A(m x n)
/** Return pair R upper triangel matrix and Q= orthogonal matrix. R and Q are always dense2D **/
template <typename Matrix, typename MatrixQ, typename MatrixR>
void qr(const Matrix& A, MatrixQ& Q, MatrixR& R)
{
    vampir_trace<4013> tracer;
    typedef typename Collection<Matrix>::value_type   		    value_type;
    typedef typename Collection<Matrix>::size_type    		    size_type;
    typedef typename Magnitude<value_type>::type      		    magnitude_type;
    typedef mtl::dense_vector<value_type, vec::parameters<> >       vector_type;
    
    size_type        ncols = num_cols(A), nrows = num_rows(A), 
                     mini= ncols == nrows ? ncols - 1 : (nrows >= ncols ? ncols : nrows);
    magnitude_type   factor= magnitude_type(2);

    Q= 1;
    for (size_type i = 0; i < mini; i++) {
	irange r(i, imax); // Intervals [i, n-1]
	vector_type   w(R[r][i]), v(householder_s(w)); 

	// R-= 2*v*(v'*R)
	MatrixR Rsub(R[r][r]);
	vector_type tmp(-factor * trans(Rsub) * v);
	rank_one_update(Rsub, v, tmp);
	
	//update Q: Q-= 2*(v*Q)*v'
	MatrixQ Qsub(Q[iall][r]);
	vector_type qtmp(-factor * Qsub * v);
	rank_one_update(Qsub, qtmp, v);
    } //end for
}

/// QR-Factorization of matrix A(m x n)
template <typename Matrix>
std::pair<mtl::mat::dense2D<typename Collection<Matrix>::value_type, mat::parameters<> >,
 	  mtl::mat::dense2D<typename Collection<Matrix>::value_type, mat::parameters<> > > 
inline qr(const Matrix& A)
{
    mtl::mat::dense2D<typename Collection<Matrix>::value_type, mat::parameters<> >  R(A), Q(num_rows(A),num_rows(A));
    qr(A, Q, R);
    return std::make_pair(Q,R);
}



// QR-Factorization of matrix A
// Return Q and R with A = Q*R   R upper triangle and Q othogonal
template <typename Matrix>
std::pair<typename mtl::mat::dense2D<typename Collection<Matrix>::value_type, mat::parameters<> >,
	  typename mtl::mat::dense2D<typename Collection<Matrix>::value_type, mat::parameters<> > >
inline qr_factors(const Matrix& A)
{
	vampir_trace<4014> tracer;
    using std::abs;
    typedef typename Collection<Matrix>::value_type   value_type;
    // typedef typename Magnitude<value_type>::type      magnitude_type; // to multiply with 2 not 2+0i
    typedef typename Collection<Matrix>::size_type    size_type;
    size_type        ncols = num_cols(A), nrows = num_rows(A);
    value_type       zero= math::zero(A[0][0]), one= math::one(A[0][0]);

    //evaluation of Q
    Matrix  Q(nrows, nrows), Qk(nrows, nrows), HEL(nrows, ncols), R(nrows, ncols), R_tmp(nrows, ncols);
    Q= one; R= zero; HEL= zero;

    boost::tie(Q, R_tmp)= qr(A);
    R= upper(R_tmp);
   
    return std::make_pair(Q,R);
}

}} // namespace mtl::matrix


#endif // MTL_MATRIX_QR_INCLUDE

