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
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/itl/pc/ilut.hpp>

template <typename Matrix>
double test_factorization(const Matrix& A, unsigned p, double tau)
{
    itl::pc::ilut<Matrix>  P(A, p, tau);
    //itl::pc::ilu_0<Matrix>  P(A);
    Matrix L(P.get_L()), U(P.get_U()), I(num_rows(A), num_cols(A));
    I= 1.0;
    L+= I;
    invert_diagonal(U);
    //std::cout << "L is\n" << L << '\n';
    //std::cout << "U is\n" << U << '\n';
    Matrix LU(L*U);
    //std::cout << "LU is\n" << LU << '\n';
    LU-= A;
    //std::cout << "LU-A is\n" << LU << '\n';
    double diff_norm= frobenius_norm(LU);
    std::cout << "|A-LU|_1 with p = " << p << ", tau = " << tau << " is " << diff_norm << '\n';
    return diff_norm;
}

int main()
{
#ifndef __PGI
    // For a more realistic example set sz to 1000 or larger
    const unsigned size = 4, N = size * size;

    mtl::compressed2D<double>          A(N, N);
    laplacian_setup(A, size, size);
       
    std::cout << "A is\n" << A << '\n';
    MTL_THROW_IF(test_factorization(A, 3, 0.001) > 0.24, mtl::logic_error("ILUT(3, 0.001) too bad"));

#if 0
    for(unsigned i= 2; i <= size; i++)
	for (double f= 0.5; f > 0.000001; f/= 2)
	    test_factorization(A, i, f);

    itl::pc::ic_0<matrix_type, float>  P(A);
    mtl::dense_vector<double>          x(N), y, yc(N);
    iota(x);
    yc= 1.03194,1.60198,1.45357,2.5258,3.55386,3.16,2.91787,4.09338,3.81335;
    
    std::cout << x << '\n';
    y= solve(P, x);
    std::cout << y << '\n';

    MTL_THROW_IF(two_norm(mtl::dense_vector<double>(y - yc)) > 0.001, 
		 mtl::logic_error("IC(0) doesn't yield expected result"));	
#endif

#endif
    return 0;
}
