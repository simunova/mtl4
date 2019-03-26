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

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;


double f(double) { cout << "double\n"; return 1.0; } 
complex<double> f(complex<double>) { cout << "complex\n"; return complex<double>(1.0, -1.0); }

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    cout << "\n" << name << "\n";

    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    typedef typename mtl::dense_vector<Scalar>            Vector;

    std::size_t size= num_cols(A);
    Matrix L(size, size), U(size, size);

    Scalar c= f(Scalar(1));   
    cout << "c is: " << c << "\n";

    for (std::size_t i= 0; i < size; i++)
	for(std::size_t j= 0; j < size; j++) {
	    U[i][j]= i <= j ? c * Scalar(i+j+2) : Scalar(0);
	    L[i][j]= i > j ? c * Scalar(i+j+1) : (i == j ? Scalar(1) : Scalar(0));
	}
    
    cout << "L is:\n" << L << "U is:\n" << U;
    A= L * U;

    Vector v(size);
    for (std::size_t i= 0; i < size; i++)
	v[i]= Scalar(i);

    Vector w( A*v );

    cout << "A is:\n" << A;

    Matrix PLU(A);

    mtl::dense_vector<size_type> Pv(size);
    lu(PLU, Pv);
    typename mtl::mat::traits::permutation<>::type P(permutation(Pv));
    
    cout << "Permuted A is \n" << Matrix(P * A);

    Matrix I(size, size);
    I= Scalar(1);

    Matrix PL(I + strict_lower(PLU)), PU(upper(PLU)), PA2(PL * PU);
    cout << "L [permuted] is:\n" << PL << "U [permuted] is:\n" << PU 
	 << "L * U [permuted] is:\n" << PA2
	 << "L * U is:\n" << Matrix(trans(P) * PA2);
 
    MTL_THROW_IF(one_norm(Matrix(trans(P) * PA2 - A)) > 0.1, mtl::runtime_error("Error in permuted LU factorization."));

    Matrix PUI(inv_upper(PU));
    cout << "inv(U) [permuted] is:\n" << PUI << "PUI * PU is:\n" << Matrix(PUI * PU);
    MTL_THROW_IF(one_norm(Matrix(PUI * PU - I)) > 0.1, mtl::runtime_error("Error in upper inversion."));

    Matrix PLI(inv_lower(PL));
    cout << "inv(L) [permuted] is:\n" << PLI << "PLI * PL is:\n" << Matrix(PLI * PL);
    MTL_THROW_IF(one_norm(Matrix(PLI * PL - I)) > 0.1, mtl::runtime_error("Error in lower inversion."));

    Matrix AI(PUI * PLI * P);
    cout << "inv(A) [inv(U) * inv(L) * P] is \n" << AI << "A * AI is\n" << Matrix(AI * A);
    MTL_THROW_IF(one_norm(Matrix(AI * A - I)) > 0.1, mtl::runtime_error("Error in inversion."));

    typename mtl::mat::traits::inv<Matrix>::type A_inv(inv(A));
    cout << "inv(A) is \n" << A_inv << "A * AI is\n" << Matrix(A_inv * A);
    MTL_THROW_IF(one_norm(Matrix(A_inv * A - I)) > 0.1, mtl::runtime_error("Error in inversion."));
}



int main(int, char**)
{
    using namespace mtl;
    std::size_t size= 4;
    
    dense2D<double>                                      dr(size, size);
    dense2D<complex<double> >                            dz(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);
    // compressed2D<double>                                 cr(size, size);

    test(dr, "Row-major dense");
    test(dz, "Row-major dense with complex numbers");
    test(dc, "Column-major dense");

    return 0;
}
