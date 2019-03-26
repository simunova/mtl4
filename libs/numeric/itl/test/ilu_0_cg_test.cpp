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

const int size = 3, N = size * size; 

template <typename Matrix>
void test()
{
    Matrix                             A;
    laplacian_setup(A, size, size);
       
    itl::pc::ilu_0<Matrix>             P(A);
    mtl::dense_vector<double>          x(N, 1.0), b(N);
    
    b = A * x;
    x= 0;
    
    itl::cyclic_iteration<double> iter(b, N, 1.e-6, 0.0, 1);
    cg(A, x, b, P, iter);

    MTL_THROW_IF(size == 3 && iter.iterations() > 4, mtl::runtime_error("Too many iterations in cg with ILU(0)"));
}

int main()
{
    test<mtl::compressed2D<double> >();
    //test<mtl::sparse_banded<double> >();

    return 0;
}
