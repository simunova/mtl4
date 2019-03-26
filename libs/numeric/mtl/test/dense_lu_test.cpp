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

// File: dense_lu_test.cpp
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

template <typename Matrix>
void singularity_test(const Matrix& A)
{
     try {
	Matrix B(lu_f(A));

    } catch (mtl::matrix_singular excp) {
	std::cout << "Exception for singularity successfully caught\n"; return;
    }
    throw "Singularity not detected";
}


int main(int , char**)
{
    using namespace mtl;
    dense2D<double, mat::parameters<tag::col_major> > A(5, 5);
    dense2D<double, mat::parameters<tag::col_major> > B(5, 5);
    dense2D<float, mat::parameters<tag::col_major> > bsp(4, 4);
    dense2D<float, mat::parameters<tag::col_major> > bsp1(4, 4);

    bsp(0,0)=0;
    bsp(0,1)=0;
    bsp(0,2)=1;
    bsp(0,3)=1;
    bsp(1,0)=2;
    bsp(1,1)=2;
    bsp(1,2)=2;
    bsp(1,3)=2;
    bsp(2,0)=1;
    bsp(2,1)=2;
    bsp(2,2)=2;
    bsp(2,3)=2;
    bsp(3,0)=1;
    bsp(3,1)=2;
    bsp(3,2)=3;
    bsp(3,3)=6;

    // Assign a three times the identity to A
    A= 3;
    A(2,3)=7.0;
    A(3,0)=7.0;
    A(4,1)=7.0;
    B=A;
    std::cout << "bsp is \n" << bsp << "\n";
    bsp1 = bsp;
    singularity_test(bsp1);
    // bsp1=lu_f(bsp1); // throws exception

    mtl::dense_vector<int> P;
    lu(bsp, P);
    std::cout << "LU(bsp) with pivoting \n" << bsp << "Permutation is " << P << "\n";
    std::cout << "LU_f(bsp) \n" << bsp1 << "\n";

    std::cout << "A is \n" << A << "\n";
    A=lu_f(A);
    std::cout << "LU_f(A) is \n" << A << "\n";
    lu(B, P);
    std::cout << "LU with pivoting is \n" << with_format(B, 5, 2) << "Permutation is " << P << "\n";

    return 0;
}
