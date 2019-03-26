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
#include <cmath>
#include <complex>
#include <vector>

#include <boost/numeric/mtl/mtl.hpp>


inline float f(float x) { return x; }
inline double f(double x) { return x; }

inline std::complex<double> f(std::complex<double> x) 
{ 
    return std::complex<double>(real(x), real(x)+1.0); 
}



int main(int, char**)
{

    typedef mtl::vec::fixed::dimension<5> fsize;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::row_major, fsize, true> >     rf;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::col_major, fsize, true> >     cf;
    rf= 1.0; cf= 1.0;

    std::cout << "rf is " << rf << '\n';
    std::cout << "cf is " << cf << '\n';

    return 0;
}
