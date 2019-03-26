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

#undef VERSION
#define VERSION 3

using namespace std;
using namespace mtl;

dense_vector<double> inline last_unit_vector(size_t n)
{
    dense_vector<double> v(n, 0.0);
    v[n-1]= 1;
    return v;
}

dense_vector<double> inline unit_vector(size_t k, size_t n)
{
    dense_vector<double> v(n, 0.0);
    v[k]= 1;
    return v;
}

#if VERSION == 1

dense2D<double> inline inverse_upper(dense2D<double> const& A)
{
    const size_t N= num_rows(A);
    assert(num_cols(A) == N); // Matrix must be square

    dense2D<double> Inv(N, N);

    for (size_t k= 0; k < N; ++k) {
	dense_vector<double> e_k(N);
	for (size_t i= 0; i < N; ++i)
	    if (i == k) 
		e_k[i]= 1.0;
	    else
		e_k[i]= 0.0;    

	dense_vector<double> res_k(N);
	res_k= upper_trisolve(A, e_k);

	for (size_t i= 0; i < N; ++i)
	    Inv[i][k]= res_k[i];
    }
    return Inv;
}

#elif VERSION == 2

dense2D<double> inverse_upper(dense2D<double> const& A)
{
    const size_t N= num_rows(A);
    assert(num_cols(A) == N); // Matrix must be square

    dense2D<double> Inv(N, N);

    for (size_t k= 0; k < N; ++k) {
	dense_vector<double> e_k(N);
	for (size_t i= 0; i < N; ++i)
	    e_k[i]= i == k ? 1.0 : 0.0;

	for (size_t i= 0; i < N; ++i)
	    Inv[i][k]= upper_trisolve(A, e_k)[i];
    }
    return Inv;
}

#elif VERSION == 3

dense2D<double> inverse_upper(dense2D<double> const& A)
{
    const size_t N= num_rows(A);
    assert(num_cols(A) == N); // Matrix must be square

    dense2D<double> Inv(N, N);
    Inv= 0;

    for (size_t k= 0; k < N; ++k) {
	irange r(0, k+1);
	Inv[r][k]= upper_trisolve(A[r][r], unit_vector(k, k+1));
    }
    return Inv;
}

#endif

dense2D<double> inline inverse_lower(dense2D<double> const& A)
{
    dense2D<double> T(trans(A));
    return dense2D<double>(trans(inverse_upper(T)));
}

dense2D<double> inline inverse(dense2D<double> const& A)
{
    assert(num_cols(A) == num_rows(A)); // Matrix must be square

    dense2D<double>          PLU(A);
    dense_vector<size_t>   Pv(num_rows(A));

    lu(PLU, Pv);
    dense2D<double>  I(num_rows(A), num_cols(A));
    I= 1;
    dense2D<double>  PU(upper(PLU)), 
        PL(strict_lower(PLU) + I);

    return dense2D<double>(inverse_upper(PU) * inverse_lower(PL) * permutation(Pv));
}



int main(int, char**)
{
    const unsigned size= 3;
    typedef dense2D<double>      Matrix;
    Matrix   A(size, size);
    A=  4, 1, 2, 
	1, 5, 3,
	2, 6, 9; 

    cout << "A is:\n" << A;

    Matrix LU(A);
    mtl::dense_vector<size_t> Pv(size);
    lu(LU, Pv);

    Matrix P(permutation(Pv));
    cout << "Permutation vector is " << Pv << "\nPermutation matrix is\n" << P;

    cout << "Permuted A is \n" << Matrix(P * A);
    Matrix I(size, size); I= 1;
    //Matrix I(mat::identity(size, size)), L(I + strict_lower(LU)), U(upper(LU));
    Matrix L(I + strict_lower(LU)), U(upper(LU));

    Matrix UI(inverse_upper(U));
    cout << "inverse(U) [permuted] is:\n" << UI << "UI * U is:\n" << Matrix(UI * U);
    assert(one_norm(Matrix(UI * U - I)) < 0.1);
    
    Matrix LI(inverse_lower(L));
    cout << "inverse(L) [permuted] is:\n" << LI << "LI * L is:\n" << Matrix(LI * L);
    assert(one_norm(Matrix(LI * L - I)) < 0.1);


    Matrix AI(UI * LI * P);
    cout << "inverse(A) [UI * LI * P] is \n" << AI << "A * AI is\n" << Matrix(AI * A);
    assert(one_norm(Matrix(AI * A - I)) < 0.1);
   
    Matrix A_inverse(inverse(A));
    cout << "inverse(A) is \n" << A_inverse << "A * AI is\n" << Matrix(A_inverse * A);
    assert(one_norm(Matrix(A_inverse * A - I)) < 0.1);

    Matrix A_e(inv(A));
    cout << "inv(A) is \n" << A_e << "\n";
    
    return 0;
}
