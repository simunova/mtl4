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
    mtl::dense2D<double> A(3, 3), B(3, 3), C(3, 3);
    A= 2.0; B= 3.0;
#if 0
    A= B * C * 8;
    cout << "B * C * 8 is\n" << A << '\n';

    A= B * (C * 8);
    cout << "B * (C * 8) is\n" << '\n';
#endif
    return 0;
}
