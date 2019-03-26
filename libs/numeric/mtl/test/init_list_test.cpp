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
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  

#ifdef MTL_WITH_INITLIST
template <typename Vector>
void test(const char* name)
{
    const Vector v= {3, 4, 5};

    mtl::io::tout << name << ", v: " << v << "\n"; 
    MTL_THROW_IF(v[0] != 3.0, mtl::runtime_error("wrong"));

    Vector w;
    w= {2, 4, 7};
    mtl::io::tout << "w: " << w << "\n"; 
    MTL_THROW_IF(w[0] != 2.0, mtl::runtime_error("wrong"));
}
#else
template <typename Vector>
void test(const char* ) {}
#endif 

int main(int , char**)
{
    using mtl::vec::parameters;
    using namespace mtl;

    test<dense_vector<float> >("test float");
    test<dense_vector<double> >("test double");
    test<dense_vector<float, parameters<row_major> > >("test float in row vector");    

    return 0;
}
 














