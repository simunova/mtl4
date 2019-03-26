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
    

template <typename Vector>
void test(Vector&, const char* name)
{
    using mtl::sum; using mtl::product; using mtl::one_norm;

    Vector a(3), b(3), res(3);
    a= 1, 2, 3; b= 4, 5, 6; res= -3, 6, -3;
    
    std::cout << name << ": cross(" << a << ", " << b << ") is " << cross(a, b) << '\n';
    
    MTL_THROW_IF(one_norm(Vector(cross(a, b) - res)) > 0.0001, mtl::runtime_error("Wrong cross product with dimension 3!\n"));

    Vector c(7), d(7);
    c= 1, 2, 3, 4, 5, 6, 7;
    d= 9, 8, 7, 6, 5, 4, 3;

    std::cout << "cross(" << c << ", " << d << ") is " << cross(c, d) << '\n';

    // What is the result???
}
 

int main(int, char**)
{
    using namespace mtl;
    using mtl::vec::parameters;

    dense_vector<float>   u;
    dense_vector<double>  x;
    dense_vector<std::complex<double> >  xc;

    std::cout << "Testing vector operations\n";

    test(u, "test float");
    test(x, "test double");

    test(xc, "test complex<double>");

    dense_vector<float, parameters<row_major> >   ur(5);
    //test(ur, "test float in row vector");

    return 0;
}
 














