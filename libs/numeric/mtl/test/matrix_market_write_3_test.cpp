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
 

using namespace std;  

std::string program_dir; // Ugly global variable !!!

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    std::cout << "\n" << name << "\n";
    A= mtl::mat::hilbert_matrix<>(4, 3);

    string fname( mtl::io::join( program_dir, string("matrix_market/write_test_3_") + string(name) + string(".mtx") ) );
    cout << "File name is " << fname << "\nA is\n" << A;
    
    mtl::io::matrix_market_ostream oms(fname);
    oms << A;
    oms.close();

    Matrix B, C; 
    B= mtl::io::matrix_market(fname);
    cout << "\nRead back results in\n" << B;

    C= A - B;
    cout << "\nDifference is\n" << C << "\none_norm of it == " << one_norm(C) << '\n';
}


int main(int, char* argv[])
{
    using namespace mtl;

    compressed2D<double>                             cdr(4, 3);
    compressed2D<float>                              cfr(4, 3);
    // compressed2D<int>                                cir(4, 3); // a Hilbert matrix as int is nonsense
    compressed2D<std::complex<double> >              ccr(4, 3);

    program_dir= mtl::io::directory_name(argv[0]);
    test(cdr, "compressed2D_double");
    test(cfr, "compressed2D_float");
    test(ccr, "compressed2D_complex");

    return 0;
}
