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
#include <vector>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>

using namespace mtl;
using namespace mtl::recursion; 
using namespace std;  


int main(int argc, char* argv[])
{
    // typedef dense2D<double> matrix_t;
    typedef morton_dense<double,  doppled_32_row_mask> matrix_t;
    const int                                          n= 10;
    matrix_t                                           A1(n, n), C1(n, n);

    // fill matrix, e.g.
    hessian_setup(A1, 1.0);
    hessian_setup(C1, 2.0);

    typedef mat::recursator<matrix_t>               recursator_t;
    recursator_t                                      ra1(A1), rc1(C1);
    std::vector<recursator_t>                         v;
    v.push_back(ra1); v.push_back(rc1);

    // Of course you can also add sub-matrices
    v.push_back(rc1.south_west());

    // Only correct with latest revision (due to subtle problem just fixed)
    for (int i= 0; i < v.size(); i++)
	cout << "Recursator " << i << ": \n" << *v[i].north_west() << "\n";

    // Should also work with older version
    for (int i= 0; i < v.size(); i++)
	cout << "Recursator " << i << ": \n", print_matrix(*v[i].north_west()), cout << "\n";

    return 0;
}
