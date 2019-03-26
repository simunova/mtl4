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

int main(int, char**)
{
    using namespace std;
    using mtl::dense2D; using mtl::dense_vector;
	
    dense2D<double> *A;
    A = new dense2D<double>[2];
    A[0] = dense2D<double>(3,3);

    A[1] = dense2D<double>(4,4);
    A[0]= 0;
    A[1]= 0;
    cout << A[0] << endl;
    cout << A[1] << endl;
    delete[] A;

    dense_vector<double> *x;
    x = new dense_vector<double>[2];
    x[0] = dense_vector<double>(3);
    x[1] = dense_vector<double>(4);
    x[0] = 0;
    x[1] = 0;
    // x[1] = dense_vector<double>(5); // Must throw an exception because 
    cout << x[0] << endl;
    cout << x[1] << endl;
    delete[] x;

    return 0;
}
