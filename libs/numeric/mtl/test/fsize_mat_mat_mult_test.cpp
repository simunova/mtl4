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


int main(int , char**)
{
    using namespace mtl;
 
    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<2, 2> > fmat_para;

    dense2D<double, fmat_para>        A, B, B2; 

    A= 1., 2.,
       3., 4.;
    B= A * A;
    std::cout << "A * A is\n" << B << "\n";
    
    B2= 7, 10,
       15, 22;
    B2-= B;
    MTL_THROW_IF(one_norm(B2) > 0.001, runtime_error("Matrix product wrong"));

    return 0;
}
