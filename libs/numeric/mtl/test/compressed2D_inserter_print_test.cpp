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

template <typename Matrix>
void test(const char* name, Matrix& A)
{
    std::cout << name << "\n";

    mtl::mat::inserter<Matrix> ins(A, 3);
    ins[0][1] << 4.; ins[0][2] << 7;
    // ins[0][3] << 6.; ins[0][0] << 8;

    ins.print(0);
    ins.print(1);
    // ins.make_empty(0);
    ins.print();
}

 
int main(int, char**)
{
    mtl::compressed2D<double> A(4, 5);
    // mtl::compressed2D<double> A(4, 5);

    test("CRS", A);


    return 0;
}
