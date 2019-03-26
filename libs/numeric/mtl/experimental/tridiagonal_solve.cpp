// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

using namespace std;

int main(int , char** ) 
{
    typedef mtl::compressed2D<double> matrix_type;
    typedef mtl::dense_vector<double> vector_type;

    const int size= 10;
    matrix_type A(size, size);

    // Set up a non-singular tridiagonal matrix
    {
	mtl::mat::inserter<matrix_type> ins(A, 3);
	for (int i= 0; i < size; i++) {
	    if (i > 0) ins[i][i-1] << -0.8;
	    ins[i][i] << 3;
	    if (i+2 < size) ins[i][i+1] << -0.8;
	}
    }
    cout << "A is\n" << A;

    vector_type x(size), b(size, 1.0);
    
    itl::pc::ilu_0<matrix_type>   P(A);
    x= solve(P, b);

    cout << "x is " << x << '\n'
	 << "A*x is " << vector_type(A*x) << '\n';

    return 0;
}
