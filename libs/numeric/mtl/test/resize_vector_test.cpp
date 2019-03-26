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

#include <boost/numeric/mtl/mtl.hpp>



using namespace std;  
 

int main(int, char**)
{
    mtl::dense_vector<float>    x;

    MTL_THROW_IF(size(x) != 0, mtl::runtime_error("vector should be empty"));

    x.change_dim(5);
    MTL_THROW_IF(size(x) != 5, mtl::runtime_error("vector should have size 5"));

    x= 3.0;
    cout << "Vector x is initialized to: " << x;

    
    x.change_dim(7);
    MTL_THROW_IF(size(x) != 7, mtl::runtime_error("vector should have size 7"));

    x= 3.0;
    cout << "Vector x after resizing (larger): " << x;

    x.change_dim(4);
    MTL_THROW_IF(size(x) != 4, mtl::runtime_error("vector should have size 4"));

    x= 3.0;
    cout << "Vector x after resizing (smaller): " << x;

    return 0;
}
 














