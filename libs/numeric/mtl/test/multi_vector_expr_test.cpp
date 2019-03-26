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

#include <cmath>
#include <iostream>
#include <vector>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;

template <typename Vector>
void test(const char* name)
{
    cout << "Testing multi_vector with " << name << endl;
    mtl::multi_vector<Vector> A(4, 6), B(4, 6);
    A= 3.0;
    cout << "A is\n" << A << endl;

    B= A;
    cout << "B= A yields\n" << B << endl;
    cout << "A is\n" << B << endl;

    B= A + A;
    cout << "A + A is\n" << B << endl;
    MTL_THROW_IF(B[1][1] != 6.0, mtl::runtime_error("Wrong value on diagonal\n"));
    MTL_THROW_IF(B[1][0] != 0.0, mtl::runtime_error("Wrong value off diagonal\n"));

    B= 2 * A;
    cout << "2 * A is\n" << B << endl;
    MTL_THROW_IF(B[1][1] != 6.0, mtl::runtime_error("Wrong value on diagonal\n"));
    MTL_THROW_IF(B[1][0] != 0.0, mtl::runtime_error("Wrong value off diagonal\n"));

    B= 2 * A + A;
    cout << "2 * A is\n" << B << endl;
    MTL_THROW_IF(B[1][1] != 9.0, mtl::runtime_error("Wrong value on diagonal\n"));
    MTL_THROW_IF(B[1][0] != 0.0, mtl::runtime_error("Wrong value off diagonal\n"));

    B= 2.0 * A + 3 * A;
    cout << "2 * A + 3 * A is\n" << B << endl;
    MTL_THROW_IF(B[1][1] != 15.0, mtl::runtime_error("Wrong value on diagonal\n"));
    MTL_THROW_IF(B[1][0] != 0.0, mtl::runtime_error("Wrong value off diagonal\n"));

    B= 2.0 * A + 3 * A - 2.6 * A;
    cout << "2 * A + 3 * A - 2.6 * A is\n" << B << endl;
    MTL_THROW_IF(std::abs(B[1][1] - 7.2) > 0.001, mtl::runtime_error("Wrong value on diagonal\n"));
    MTL_THROW_IF(B[1][0] != 0.0, mtl::runtime_error("Wrong value off diagonal\n"));

    mtl::multi_vector<Vector> C;
    C.change_dim(3, 10);
    C= 7;
    cout << "C is\n" << C;
}

int main(int, char**)
{
    test<mtl::dense_vector<double> >("dense_vector<double>");

    return 0;
}
