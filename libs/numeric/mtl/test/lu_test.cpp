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
void singularity_test1(const Matrix& A)
{
    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    try {
	Matrix B(A);
	B[mtl::iall][0]= Scalar(0);
	// cout << "B is:\n" << B;
	lu(B);
    } catch (mtl::matrix_singular excp) {
	cout << "Exception 1 for singularity successfully caught\n"; return;
    }
    throw "Singularity (test1) not detected";
}

template <typename Matrix>
void singularity_test2(const Matrix& A)
{
    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    try {
	Matrix B(A);
	B[mtl::iall][0]= Scalar(0);
	mtl::dense_vector<int> p;
	lu(B, p);
	cout << "B is:\n" << B << endl;
    } catch (mtl::matrix_singular excp) {
	cout << "Exception 2 for singularity successfully caught\n"; return;
    }
    throw "Singularity (test2) not detected";
}

template <typename Matrix>
void singularity_test3(const Matrix& A)
{
    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    try {
	Matrix B(A);
	B[mtl::iall][0]= Scalar(0);
	B[0][0]= Scalar(1e-15);
	lu(B, 2e-15);
    } catch (mtl::matrix_singular excp) {
	cout << "Exception 3 for singularity successfully caught\n"; return;
    }
    throw "Singularity (test3) not detected";
}

template <typename Matrix>
void singularity_test4(const Matrix& A)
{
    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    try {
	Matrix B(A);
	B[mtl::iall][0]= Scalar(0);
	B[3][0]= Scalar(1e-15);
	mtl::dense_vector<int> p;
	lu(B, p, 2e-15);
    } catch (mtl::matrix_singular excp) {
	cout << "Exception 4 for singularity successfully caught\n"; return;
    }
    throw "Singularity (test4) not detected";
}




template <typename Matrix>
void test(Matrix& A, const char* name)
{
    cout << "\n" << name << "\n";

    typedef typename mtl::Collection<Matrix>::value_type  Scalar;
    typedef typename mtl::dense_vector<Scalar>            Vector;

    unsigned size= unsigned(num_cols(A));
    Matrix L(size, size), U(size, size);

    Scalar c= f(Scalar(1));   
    cout << "c is: " << c << "\n";

    for (unsigned i= 0; i < size; i++)
	for(unsigned j= 0; j < size; j++) {
	    U[i][j]= i <= j ? c * Scalar(i+j+2) : Scalar(0);
	    L[i][j]= i > j ? c * Scalar(i+j+1) : (i == j ? Scalar(1) : Scalar(0));
	}
    
    cout << "L is:\n" << L << "U is:\n" << U;
    A= L * U;

    Vector v(size);
    for (unsigned i= 0; i < size; i++)
	v[i]= Scalar(i);

    Vector w( A*v );

    cout << "A is:\n" << A;

    Matrix LU(A);
    lu(LU);
    cout << "LU decomposition of A is:\n" << LU;

    Matrix I(size, size);
    I= Scalar(1);
    Matrix tmp(I + strict_lower(LU)), A2(tmp * upper(LU));
    cout << "L * U is:\n" << A2;

    Matrix B( lu_f(A) );

    Vector v2( upper_trisolve(upper(LU), unit_lower_trisolve(strict_lower(LU), w)) );

    cout << "LU decomposition of A (as function result) is:\n" << B;
    cout << "upper(LU) is:\n" << upper(LU) << "strict_lower(LU) is:\n" << strict_lower(LU);
    cout << "v2 is " << v2 << "\n";

    MTL_THROW_IF(abs(v[1] - v2[1]) > 0.1, mtl::runtime_error("Error using tri_solve"));

    Vector v3( lu_solve_straight(A, w) );
    MTL_THROW_IF(abs(v[1] - v3[1]) > 0.1, mtl::runtime_error("Error in solve"));

    Vector v4( lu_solve(A, w) );
    cout << "v4 is " << v4 << "\n";
    MTL_THROW_IF(abs(v[1] - v4[1]) > 0.1, mtl::runtime_error("Error in solve"));

    singularity_test1(A);
    singularity_test2(A);
    singularity_test3(A);
    singularity_test4(A);

    mtl::mat::lu_solver<Matrix> lus(A);
    Vector v5(size);
     
    lus.solve(w, v5);
    cout << "v5 is " << v5 << "\n";
    MTL_THROW_IF(abs(v[1] - v5[1]) > 0.1, mtl::runtime_error("Error in solve"));
}



int main(int, char**)
{
    using namespace mtl;
    unsigned size= 4;
    
    dense2D<double>                                      dr(size, size);
    dense2D<complex<double> >                            dz(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);

    test(dr, "Row-major dense");
    test(dz, "Row-major dense with complex numbers");
    test(dc, "Column-major dense");

    return 0;
}
