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
#include <boost/numeric/itl/itl.hpp>


int main()
{
    using mtl::srange; using mtl::imax;
    // For a more realistic example set sz to 1000 or larger
    const int size = 3, N = size * size; 

    typedef mtl::compressed2D<double>  	   matrix_type;
    typedef itl::pc::ic_0<matrix_type> 	   ic_type;
    typedef itl::pc::diagonal<matrix_type> dia_type;

    mtl::compressed2D<double>          A;
    laplacian_setup(A, size, size);
    // mtl::io::tout << "A is\n" << A << '\n';

    itl::pc::concat<ic_type, dia_type, matrix_type> L(A);
    // ic_type                                         L(A);
    itl::pc::identity<matrix_type>                  R(A);

    mtl::dense_vector<double>          x(N, 1.0), b(N);
    
    b = A * x;
    x= 0;

    itl::cyclic_iteration<double> iter(b, N, 1.e-6, 0.0, 3);
    cg(A, x, b, L, iter);

    // Test if adjoint works (do not try bicg, it convergences poorly)
    x= 0;
    itl::cyclic_iteration<double> iter2(b, N, 1.e-6, 0.0, 5);
    qmr(A, x, b, L, R, iter2);

    return 0;
}
