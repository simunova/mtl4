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
#include <utility>
#include <cmath>
// #include <boost/test/minimal.hpp>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/itl/smoother/gauss_seidel.hpp>

using namespace std;  
   
int main(int, char**)
{
    using namespace mtl;

    const int s= 10;
    typedef mtl::dense_vector<double> Vector;
    typedef mtl::compressed2D<double> Matrix;
    Vector       x(s*s, 8), b(s*s);
    Matrix   A;
    laplacian_setup(A, s, s);
    if (s < 10)
      std::cout<< "x= " << x << "\n";
    
    b= A*x;
    x= 0;

    itl::gauss_seidel<Matrix> gs(A);
    for (int i =0 ; i< 30; i++)
        gs(x, b);
    
    if (s < 10) {
      std::cout<< "x=" << x << "\n";
      Vector tmp(b-A*x);
      assert(two_norm(tmp) < 1.0e-4);
    }
    
    return 0;
}
 














