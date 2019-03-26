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



int main()
{
    using namespace mtl;
    typedef dense2D<double>  mat;
    typedef dense_vector<double, mtl::parameters<mtl::tag::row_major> >  row_vector;
    typedef operations::bracket_proxy<mat, mat&, double&> proxy_type;

    mat    A(2, 3);
    A= 1, 2, 3,
       4, 5, 6;
	
    row_vector   v(A[1][iall]);
    v[1]= 7;

    std::cout << "A is\n" << A;

    proxy_type row_1(A[1]);
    v[1]= 8;

    std::cout << "A is\n" << A;
    
    return 0;
}
