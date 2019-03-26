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


using namespace std;  
    

template <typename Vector>
void test(Vector& v, const char* name)
{
    typedef typename mtl::Collection<Vector>::value_type value_type;
    using mtl::max;

    v = value_type(-2);
    std::cout << "\n" << name << "  --- v = " << v << std::endl;
    std::cout << "max(v) is " << max(v) << std::endl;

    if (max(v) >= value_type(-1))
	throw "Max value too large\n";

}
 

int main(int, char**)
{
    using namespace mtl;

    dense_vector<float>   u(5);
    dense_vector<short int>     vs(5);
    dense_vector<int>     v(5);
    dense_vector<long int>     vl(5);
    dense_vector<double>  x(5);

    test(vs, "test short int");
    test(v, "test int");
    test(vl, "test long int");
    test(u, "test float");
    test(x, "test double");

    dense_vector<float, vec::parameters<row_major> >   ur(5);
    test(ur, "test float in row vector");

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    dense_vector<long double>  xl(5);
    test(xl, "test long double");
#endif   

    return 0;
}
 














