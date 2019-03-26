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

#ifndef MTL_MATRIX_SVD_INCLUDE
#define MTL_MATRIX_SVD_INCLUDE

#include <cmath>
#include <limits>
#include <algorithm>
#include <boost/numeric/mtl/matrix/strict_upper.hpp>
#include <boost/numeric/mtl/operation/diagonal.hpp>
#include <boost/numeric/mtl/operation/one_norm.hpp>
#include <boost/numeric/mtl/operation/sub_matrix.hpp>
#include <boost/numeric/mtl/operation/trans.hpp>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

/// Returns A=S*V*D' for matrix A as references
template <typename Matrix>
inline void svd(const Matrix& A, Matrix& S, Matrix& V, Matrix& D, double tol= 10e-10)
{
	vampir_trace<3037> tracer;
    typedef typename Collection<Matrix>::value_type   value_type;
    typedef typename Collection<Matrix>::size_type    size_type;
    size_type        ncols= num_cols(A), nrows= num_rows(A), loops, col= ncols, row= nrows;
    value_type       ref, zero= math::zero(ref), one= math::one(ref); 
    double 	     err(std::numeric_limits<double>::max()), e, f;

    if (nrows != ncols) // important for right dimension
	std::swap(row, col);
    
     //init
    Matrix Q(row,row),  R(row,col),  VT(row,col), E(row,col), 
 	   QT(col,col), RT(col,row);

    loops= 100 * std::max(nrows,ncols);
    S= one; D= one; E= zero;
    for (size_type i= 0; err > tol && i < loops; ++i) {
	boost::tie(QT, RT)= qr(V);
 	S*= QT;
	VT= trans(RT);
	boost::tie(Q, R)= qr(VT);
	D*= Q;
	E= triu(R,1);
	V= trans(R);

	//ready for exit when upper(R)=0
	f= two_norm(diagonal(R));
	e= one_norm(E);
	if ( f== zero ) f= 1;
	err= e/f;
    } //end for
    
    {
	V= 0;  
	mtl::mat::inserter<Matrix>  ins_V(V);
	mtl::mat::inserter<Matrix,  mtl::operations::update_times<value_type> > ins_S(S);

	for (size_type i= 0, end= std::min(nrows, ncols); i < end; i++) {
	    ins_V[i][i] << std::abs(R[i][i]);
	    if (R[i][i] < zero) 	
		for (size_type j= 0; j < nrows; j++) 
		    ins_S[j][i] << -1;  //carefull changing: multiplication with minus one
	}
    }
}

/// Returns A=S*V*D' for matrix A as triplet
template <typename Matrix>
boost::tuple<Matrix, Matrix, Matrix >
inline svd(const Matrix& A, double tol= 10e-10)
{
	vampir_trace<3038> tracer;
    typedef typename Collection<Matrix>::size_type    size_type;
    size_type    ncols= num_cols(A), nrows= num_rows(A), col= ncols, row= nrows;
    if (nrows != ncols) // important for right dimension
	std::swap(row, col);

    Matrix       ST(col,col), V(A), D(row,row);
    svd(A, ST, V, D, tol);
    return boost::make_tuple(ST, V, D);
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_SVD_INCLUDE
