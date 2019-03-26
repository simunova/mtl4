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


int main(int, char**)
{
    typedef mtl::dense2D<double> Matrix;

    Matrix A(3,4);
    A[0][0]=1;A[0][1]=2;A[0][2]=3;A[0][3]=4;
    A[1][0]=5;A[1][1]=6;A[1][2]=7;A[1][3]=8;
    A[2][0]=9;A[2][1]=10;A[2][2]=11;A[2][3]=12;

    std::cout << "First Matrix A is: \n" << A << "\n";

    Matrix B = mtl::clone(sub_matrix(A,0,1,0,1));
    std::cout << "Matrix B is: \n" << B << "\n";

    Matrix C = mtl::clone(sub_matrix(A,0,1,0,2));
    std::cout << "Matrix C is: \n" << C << "\n";

    Matrix D = mtl::clone(sub_matrix(A,0,1,0,3));
    std::cout << "Matrix D is: \n" << D << "\n";

    Matrix E = mtl::clone(sub_matrix(A,0,1,0,4));
    std::cout << "Matrix E is: \n" << E << "\n";

    Matrix f = mtl::clone(sub_matrix(A,0,2,0,1));
    std::cout << "Matrix f is: \n" << f << "\n";

    MTL_THROW_IF(f[1][0] != 5.0, mtl::runtime_error("Clown has wrong value."));
    f[1][0]= 6.0;
    std::cout << "Matrix f after f[1][0]= 6.0 is: \n" << f << "\n";
    std::cout << "Matrix A after f[1][0]= 6.0 is: \n" << A << "\n";

    MTL_THROW_IF(f[1][0] != 6.0, mtl::runtime_error("Clown has wrong value after change."));
    MTL_THROW_IF(A[1][0] != 5.0, mtl::runtime_error("Original matrix was changed."));
    

    Matrix g = mtl::clone(sub_matrix(A,0,3,0,1));
    std::cout << "Matrix g is: \n" << g << "\n";

    Matrix h = mtl::clone(sub_matrix(A,0,3,0,1));
    std::cout << "Matrix h is: \n" << h << "\n";

    Matrix i = mtl::clone(sub_matrix(A,0,2,0,4));
    std::cout << "Matrix i is: \n" << i << "\n";

    Matrix k = mtl::clone(sub_matrix(A,0,2,0,2));
    std::cout << "Matrix k is: \n" << k << "\n";

    Matrix l = mtl::clone(sub_matrix(A,1,3,0,2));
    std::cout << "Matrix l is: \n" << l << "\n";
  

    return 0;
}
