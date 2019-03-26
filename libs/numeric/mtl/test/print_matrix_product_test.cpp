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
    
    double aa[2][2]= {{1., 2.},
		      {3., 4.}};
    mtl::dense2D<double> A(aa), B(2,2);
    B = A*A;
    
    cout << (A*A) << endl << B << endl;
    
    MTL_THROW_IF(B(0, 1) != (A*A)(0,1), mtl::runtime_error("Wrong value in matrix product expression!\n"));

    return 0;
}
