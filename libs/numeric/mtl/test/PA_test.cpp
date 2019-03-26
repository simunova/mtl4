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
    // Test that matrix product don't throw index ouf of range (anymore)
    using namespace std;
    
    double A_array[][1]= {{1.}, {2.,}};
    mtl::dense2D<double, mtl::mat::parameters<> >  A(A_array);
    std::cout << "A=\n" << A << std::endl;

    mtl::dense_vector<unsigned, mtl::vec::parameters<> >  perm_v(2);
    perm_v= 1, 0;
 
    mtl::dense2D<double, mtl::mat::parameters<> > PA(permutation(perm_v)*A);
    std::cout << "PA=\n" << PA << std::endl;

    return 0;
}
