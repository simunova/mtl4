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
    mtl::dense2D<double> A(mtl::mat::hilbert_matrix<>(5, 6));
    A[2][0]= 1234.5; A[2][2]= 0.0004; A[2][4]= 1234567;

    cout << "Pretty print test. The following matrix should be printed prettily with aligned columns.\n" << A;
    
    return 0;
}
 
