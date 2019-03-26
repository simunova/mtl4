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

#ifndef __PGI
#include <boost/lambda/lambda.hpp>
#endif

int main(int , char**)
{
#ifndef __PGI
    using namespace std;
    using boost::lambda::_1;

    mtl::compressed2D<double> A;
    laplacian_setup(A, 3, 4);

    assign_each_nonzero(A, constant_ref(_1) * constant_ref(_1));
    mtl::io::tout << "A is\n" << A << std::endl;
    
    MTL_THROW_IF(abs(A[0][0] - 16.0) > 0.001, mtl::runtime_error("Wrong value in diagonal"));
    MTL_THROW_IF(abs(A[0][1] - 1.0) > 0.001, mtl::runtime_error("Wrong value in off-diagonal"));
    MTL_THROW_IF(abs(A[0][2]) > 0.001, mtl::runtime_error("Entry should be empty"));
#endif

    return 0;
}
