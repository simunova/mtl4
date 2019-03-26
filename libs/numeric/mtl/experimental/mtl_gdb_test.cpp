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
#include <vector>
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;
using namespace mtl;

int main(int, char**) 
{
    float A[2][3]= {{2, 3, 6}, {1, 2, 0}};
    A[0][0]= 24.;

    std::complex<float> z(2, 5);

    std::vector<std::complex<float> > v(3, z);

    std::vector<float> v2(3, 2.5);

    dense_vector<double>   w(3), w2(9);
    w= 4, 5, 6;
    iota(w2);

    dense2D<double>  B(2, 3);
    B= 2;

    dense2D<double, mat::parameters<col_major> >  C(2, 3);
    C= 2, 3, 4,
       5, 6, 7;

    cout << "C is\n" << C;

    compressed2D<double>  D(2, 3);
    D= 2;

    compressed2D<double, mat::parameters<col_major> >  E(2, 3);
    E= 2, 3, 4,
       5, 6, 7;

    cout << "E is\n" << E;
     
    compressed2D<float> F;
    laplacian_setup(F, 3, 4);    
    cout << "F is\n" << F;

    return 0;
}
