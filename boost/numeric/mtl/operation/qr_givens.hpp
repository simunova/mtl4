// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

// Author: Marc Hartung

#ifndef MTL_MATRIX_QR_GIVENS_INCLUDE
#define MTL_MATRIX_QR_GIVENS_INCLUDE

#include <cmath>
#include <utility>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>


namespace mtl { namespace mat {


/// General Given's transformator
/**
 * Input-matrix A will be splitted in A = Q*R. Q and R are accessible from getQ() and getR().
 * \sa givens
 */    
template <typename Matrix>
class qr_givens_solver {
    
    typedef typename Collection<Matrix>::value_type value_type;
    typedef typename Collection<Matrix>::size_type size_type;

public:
    
   /** \brief Constructor needs a sqare-matrix as input
     * \param MTL-Matrix type
     *  
     */
    
    qr_givens_solver(const Matrix& IN) : R(IN), G(2,2), Q(num_cols(IN), num_rows(IN)) {
	value_type one = math::one(c);
	Q = one;
	eps = 1.0e-8;
    }
    
    /** \brief Returns Q as the input-matrix-type 
     * 
     * \return Matrix Q
     */
    
    Matrix& getQ() {
	return Q;
    }
    
    /** \brief Returns R as the input-matrix-type 
     * 
     * \return Matrix R
     */    
    Matrix& getR() {
	return R;
    }
    
    /** \brief Sets the distance to zero
     * 
     * \param value_type tolerance to zero
     */    
    void setTolerance(value_type tol) {
	eps = tol;
    }
    
    /** \brief Starts QR decompostion
     */	
    void calc() 
    {
	value_type zero= math::zero(c);
	Matrix RIter(2, num_rows(R)), Iter(2, num_rows(R));
	for(size_type i= 0; i < R.num_cols() - 1; i++) {
	    irange r(i, num_cols(R));
 	    for(size_type j= num_rows(R)-1;j>i;j--) {
		if (std::abs(R[j][i]) > eps) {
		    set_rotation(R[i][i], R[j][i]);
		    Iter[0][r] = R[i][r];
		    Iter[1][r] = R[j][r];
		    RIter[iall][r] = G*Iter[iall][r];
		    R[i][r] = RIter[0][r];
		    R[j][r] = RIter[1][r];
		    Iter[0][iall] = Q[i][iall];
		    Iter[1][iall] = Q[j][iall];
		    RIter = G*Iter;
		    Q[i][iall] = RIter[0][iall];
		    Q[j][iall] = RIter[1][iall];
		}
		R[j][i] = zero;
	    }
	    
	}
	
	for(size_type i= 0; i < num_cols(R) - 1; i++) 
	    R[i+1][i] = zero;
	
    }
        
  private:
    
    Matrix     R, G, Q;
    value_type c, s, eps;
    
    template <typename T>
    inline T square(T x) const { return x * x; }

    /// Coefficients for Given's rotation    
    void set_rotation(value_type a, value_type b)
    {
	using std::abs;
	using mtl::conj;
	
	value_type zero= math::zero(a), one= math::one(b), t;
	
	if ( b == zero ) {
	    c= one; s= zero;
	} else if ( abs(b) > abs(a) ) {
	    t = a / abs(b);
	    s = one/ sqrt(one + square(abs(t)));
	    c = s * t;
	    s *= (b/abs(b));
	} else {
	    t = b / abs(a);
	    c = one / sqrt(one + square(abs(t)));
	    s = c * t;
	    c = (a/abs(a))*c;
	}
	
	G= conj(c), conj(s),
	  -s, c;
    }
    
};

/// QR-Factorization of matrix A(m x n) based on Givens' rotation
template <typename Matrix>
std::pair<Matrix, Matrix>
inline qr_givens(const Matrix& A)
{
    qr_givens_solver<Matrix> solver(A);
    solver.calc();
    return std::make_pair(solver.getQ(), solver.getR());
}

}} // namespace mtl::matrix

#endif // MTL_MATRIX_QR_GIVENS_INCLUDE
