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
//
// Author Marc Hartung

#ifndef MTL_MATRIX_EIGENVALUE_INCLUDE
#define MTL_MATRIX_EIGENVALUE_INCLUDE

#include <cmath>
#include <complex>
#include <utility>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/qr_givens.hpp>
#include <boost/numeric/mtl/operation/misc.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>

namespace mtl { namespace mat {


/// Solver class for general eigenvalues
/** Not yet tested for complex matrices. **/
template <typename Matrix>
class eigenvalue_solver {
    
    typedef typename Collection<Matrix>::value_type   value_type;
    typedef typename Collection<Matrix>::size_type    size_type;  
        
public:
    
    /** \brief Constructor needs a sqare-matrix as input
     * \param MTL-Matrix type
     *  
     */
    eigenvalue_solver(const Matrix& IN) : ncols(num_cols(IN)), nrows(num_rows(IN)) {
	zero= math::zero(IN[0][0]);
	one= math::one(IN[0][0]);
	R = hessenberg(IN);
	I = Matrix(ncols,nrows);
	I = one;
	
	//Standartwert Initialisierung:
	
	maxIteration = 20*ncols*ncols;
	eps = 1.0e-8;
    }
    
    
    /** \brief Changes the zero-tolerance
     *  \param value_type of allowed distance to zero
     */
    void setTolerance(value_type in) {
	eps = in;
    }
    
    /** \brief Changes the number of allowed iterations
     * \param size_type variable of maximum allowed iterations
     */    
    void setMaxIteration(size_type in) {
	maxIteration = in;
    }
    
    /** \brief Starts the calculation of eigenvalues.
     */
    void calc() {
	value_type singleCont;
	std::complex<value_type> compCont;
	size_type size, allIt = 0;
	for(size_type i=ncols-1;i>0 && i<ncols && allIt<=maxIteration;i--) {
	    size = i+1; // verkleinert die Matrix-Dimensionen
	    irange r(0,size);
	    singleCont = std::abs(R[i][i-1])+1;
	    compCont = getComplex2x2EW(i) + 1;
	    while( allIt < maxIteration  && std::abs(R[i][i-1])>eps) {
		if(isRealEW(i)) { 
		    if(std::abs(R[i][i-1])<singleCont) { //SingleShift
			singleCont = std::abs(R[i][i-1]);
			singleShift(r);
		    }
		    else { //DoubleShift
			singleCont = std::abs(R[i][i-1]);
			doubleShift(r);
		    }
		    allIt++;
		}
		else {
		    if(std::abs(compCont-getComplex2x2EW(i))>eps) {
			compCont = getComplex2x2EW(i);
			doubleShift(r);
			allIt++;
		    }
		    else {
			i--;
			break;
		    }		    
		}
 	    }
	}
    }
    
    
    /** \brief Returner for the calculated eigenvalues
     * 
     * Before using get_eigenvalues() you have to use calc()!
     * \param dense_vector<complex<value_type>> of eigenvalues
     */    
    dense_vector<std::complex<value_type> > get_eigenvalues() 
    {
	using mtl::conj;
	dense_vector<std::complex<double> > res(ncols, 0.0);
	size_type i;
	for(i=ncols-1;i>0 && i<ncols;i--) {
	    if(std::abs(R[i][i-1])<eps || isRealEW(i)) { // wenn ein reeller EW
		res[i] = R[i][i];
	    }
	    else { // wenn zwei konjungiert komplexe EW

		std::complex<value_type> ews = getComplex2x2EW(i);
		res[i-1] = ews;
		res[i] = conj(ews);
		i--;
	    }
	}
	if(i==0) { // zur Vermeidung von Zugriffsfehler
	    res[0] = R[0][0];
	}    
	return res;
    }
          
    
    
private:
    Matrix R, I;
    size_type ncols, nrows, maxIteration;
    value_type zero, one, eps;
    
    
    bool isRealEW(size_type k) {
	if(std::abs(sqrt(square(R[k-1][k-1]+R[k][k])/4.0+R[k-1][k]*R[k][k-1]-R[k-1][k-1]*R[k][k])) >= 0.0) {
	    return true;
	}
	return false;
    }
    
    /** \brief Calculates real eigenvalues of a 2x2-matrix defined by col/row k arround the diagonal of the input-matrix
     * 
     * First entry of the pair is the eigenvalue closer to the (kxk)-entry (Wilkinson-shift for single shifting)
     */
    std::pair<value_type, value_type> get2x2EW(const size_type k) {
	std::pair<value_type, value_type> res;
	value_type front, back, comparator;
	
	front = (R[k-1][k-1]+R[k][k])/2.0;
	back = sqrt(square(R[k-1][k-1]+R[k][k])/4.0+R[k-1][k]*R[k][k-1]-R[k-1][k-1]*R[k][k]);
	comparator = R[k][k]-front;
	
	if( std::abs(comparator-back) < std::abs(comparator+back) ) {
	    res.first = front+back;
	    res.second = front-back;
	}
	else {
	    res.first = front-back;
	    res.second = front+back;
	}
	return res;
    }
    
    
    /** \brief Calculates a complex eigenvalue of a 2x2-matrix defined by col/row k arround the diagonal of the input-matrix
     * 
     */
    std::complex<value_type> getComplex2x2EW(const size_type k) {
	std::complex<value_type> res;
	res = square(R[k-1][k-1]+R[k][k])/4.0+R[k-1][k]*R[k][k-1]-R[k-1][k-1]*R[k][k];
	res = sqrt(res);
	res += (R[k-1][k-1]+R[k][k])/2.0;
	return res;
    }
    
    /** \brief Performes a double shift for submatrix defined by range r
     */
    
    void doubleShift(irange r) 
    {
	using mtl::conj;
	value_type s,t;
	
	if(isRealEW(r.finish()-1)) {
	    std::pair<value_type, value_type> ews = get2x2EW(r.finish()-1);
	    s = ews.first+ews.second;
	    t = ews.first*ews.second;
	}    
	else {
	    std::complex<value_type> ew = getComplex2x2EW(r.finish()-1);
	    s = 2.0*ew.real();
	    t = std::abs(ew*conj(ew));
	}
	
	Matrix RIter(R[r][r]*R[r][r] - s*R[r][r] + t*I[r][r]);
	
	qr_givens_solver<Matrix> QR(RIter);
	QR.setTolerance(eps);
	QR.calc();
	
	RIter = R[r][r];
	R[r][r] = (QR.getQ()) * RIter * trans(QR.getQ());
	
    }
    
    
    
     /** \brief Performes a single shift for submatrix defined by range r
     */
    void singleShift(irange r) {
	value_type sh = get2x2EW(r.finish()-1).first;
	if(!(std::abs(sh) >= 0.0)) {
	    sh = R[r.finish()-1][r.finish()-1];
	}
	Matrix RIter(R[r][r] - sh*I[r][r]);		//Shift wird abgezogen
	
	qr_givens_solver<Matrix> QR(RIter);		//QR-Zerlegung
	QR.setTolerance(eps);
	QR.calc();
	
	R[r][r] = ((QR.getR())*trans(QR.getQ())) + sh*I[r][r];	//QR-Iteration mit dem Rückshift
	
    }

};

/// Calculation of eigenvalues of general matrix A
/** Not yet tested for complex matrices. **/
template <typename Matrix>
dense_vector<std::complex<typename Collection<Matrix>::value_type> >
inline eigenvalues(const Matrix& A)
{
    eigenvalue_solver<Matrix> s(A);
    s.calc();
    return s.get_eigenvalues();
} 


}} // namespace mtl::matrix

#endif // MTL_MATRIX_EIGENVALUE_INCLUDE
