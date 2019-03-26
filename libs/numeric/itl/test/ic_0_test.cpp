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
#include <boost/numeric/itl/itl.hpp>

int main()
{
    // For a more realistic example set sz to 1000 or larger
    const int size = 3, N = size * size;

    typedef mtl::compressed2D<double>  matrix_type;
    mtl::compressed2D<double>          A(N, N);
    laplacian_setup(A, size, size);
       
    itl::pc::ic_0<matrix_type, float>  P(A);
    mtl::dense_vector<double>          x(N), y, yc(N);
    iota(x);
    yc= 1.03194,1.60198,1.45357,2.5258,3.55386,3.16,2.91787,4.09338,3.81335;
    
    std::cout << x << '\n';
    y= solve(P, x);
    std::cout << y << '\n';

    MTL_THROW_IF(two_norm(mtl::dense_vector<double>(y - yc)) > 0.001, 
		 mtl::logic_error("IC(0) doesn't yield expected result"));	
    return 0;
}
