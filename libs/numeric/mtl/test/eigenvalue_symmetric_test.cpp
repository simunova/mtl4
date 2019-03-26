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

// With contributions from Cornelius Steinhardt

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;


int main(int , char**)
{
    using namespace mtl;

    dense_vector<double>                    eig;

    double array[][4]= {{1,  1,   1,  0},
                        {1, -1,  -2,  0},
                        {1, -2,   1,  0},
                        {0,  0,   0, 10}};
    dense2D<double> A(array);
    std::cout << "A=\n" << A << "\n";

    eig= eigenvalue_symmetric(A,22);

    std::cout<<"eigenvalues  ="<< eig <<"\n";

    return 0;
}



