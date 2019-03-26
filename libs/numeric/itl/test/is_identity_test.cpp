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


template <typename PC>
void test(const PC& pc, bool expected)
{
    if (is_identity(pc) != expected) {
	std::cerr << "is_identity should be " << (expected ? "true" : "false") << '\n';
	throw "wrong result";
    }
}


int main()
{
    // For a more realistic example set sz to 1000 or larger
    const int size = 2, N = size * size; 

    typedef mtl::compressed2D<double>  matrix_type;
    mtl::compressed2D<double>          A(N, N), dia(N, N);
    laplacian_setup(A, size, size);

    itl::pc::identity<matrix_type>     identity(A);
    itl::pc::ic_0<matrix_type>         ic_0(A);
    itl::pc::ilu_0<matrix_type>        ilu_0(A);
    itl::pc::diagonal<matrix_type>     diagonal(A);

    test(identity, true);
    test(ic_0, false);
    test(ilu_0, false);
    test(diagonal, false);

    return 0;
}
