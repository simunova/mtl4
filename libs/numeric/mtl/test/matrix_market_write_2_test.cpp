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
#include <cmath>
#include <string>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>
 

using namespace std;  


std::string program_dir; // Ugly global variable !!!

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    laplacian_setup(A, 3, 4);
    cout << name << ":\n A is\n" << A;
    mtl::io::matrix_market_ostream oms(mtl::io::join(program_dir, "matrix_market/laplace_3x4.mtx"));
    oms << A;
}


int main(int, char* argv[])
{
    using namespace mtl;

    compressed2D<double>                             cdr;

    program_dir= mtl::io::directory_name(argv[0]);

    test(cdr, "compressed2D_double");

    return 0;
}
