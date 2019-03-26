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

#include <boost/numeric/mtl/mtl.hpp>
#include <iostream>

int main(int ,char**)
{
    using namespace mtl::mat;
    
    const unsigned n= 100;
    dense2D<double>                            A(n, n);
    morton_dense<double, mtl::doppled_64_row_mask>  C(n, n);

    A[0][0]= 71.; C[0][0]= 72.;

    double *ap= &A[0][0];
    std::cout << *ap << '\n';
    MTL_THROW_IF(*ap != 71, mtl::runtime_error("wrong value"));

    // the direct access of the first element should only be used when absolutely necessary 
    double *ap2= A.elements();
    std::cout << *ap2 << '\n';
    MTL_THROW_IF(*ap2 != 71, mtl::runtime_error("wrong value"));

    return 0;
}
