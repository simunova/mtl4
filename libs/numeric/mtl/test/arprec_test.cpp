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
    using namespace std;
    
#ifdef MTL_HAS_ARPREC

    mp::mp_init(40); 

    // mp_complex z= 0, z2= 3, z3;
    // z3= z * z2;
    // std::cout << "0 * 3 as complex is " << z3.real << "\n";

    mp_real s0 = 1.0;
    mp_real s1 = 1-s0;
    mtl::dense_vector<mp_real> v0(2,2.0);
    (1-s0)*v0;

    {
	mp_real s1 = 1-s0;
    }

    mtl::dense_vector<mp_real> v1(2,3.0);
    v1 = s1*v0;
    v1 = (1-s0)*v0;
    std::cout << "s1*v0= " << s1*v0 << "\n";
    std::cout << "(1-s0)*v0 " << (1-s0)*v0 << "\n";
    // 1-s0; memory leak !!!!!

    mp::mp_finalize(); 

#else
    cout << "ARPREC not active, test not performed.\n";
#endif

    return 0;
}
