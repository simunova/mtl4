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

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/operation/norms.hpp>
#include <boost/numeric/mtl/operation/sum.hpp>
#include <boost/numeric/mtl/operation/product.hpp>
#include <boost/numeric/mtl/operation/unroll.hpp>


using namespace std;  
    
template <typename Vector, typename Value>
void extra_test(Vector& u, Vector& v, Value)
{
    u+= 4;
    std::cout << "u+= 4 = " << u << "\n"; 

    u= v + 4;
    std::cout << "u= v + 4 = " << u << "\n"; 

    u= 4 + v;
    std::cout << "u= 4 + v = " << u << "\n"; 

    u-= 4;
    std::cout << "u-= 4 = " << u << "\n"; 

    u= v - 4;
    std::cout << "u= v - 4 = " << u << "\n"; 

    u= 4 - v;
    std::cout << "u= 4 - v = " << u << "\n"; 

    v= min(4, u);
    std::cout << "v= min(4, u) = " << v << "\n"; 

    v= min(u, 4);
    std::cout << "v= min(u, 4) = " << v << "\n"; 

    v= max(4, u);
    std::cout << "v= max(4, u) = " << v << "\n"; 

    v= max(u, 4);
    std::cout << "v= max(u, 4) = " << v << "\n"; 

    u= max(v - 4, 1);
    std::cout << "u= max(v - 4, 3) is " << u << "\n"; 
}

// not complex
template <typename Vector, typename Value>
void extra_test(Vector&, Vector&, std::complex<Value>) {}


template <typename Vector>
void test(Vector& v, const char* name)
{
    typedef typename mtl::Collection<Vector>::value_type value_type;
    typedef typename mtl::Collection<Vector>::size_type  size_type;

    using mtl::sum; using mtl::product; using mtl::one_norm; using mtl::unroll;

    for (size_type i= 0; i < size(v); i++)
	v[i]= value_type(double(i+1) * pow(-1.0, int(i))); // Amb. in MSVC

    std::cout << "\n" << name << "  --- v = " << v << std::endl;

    Vector w(v + v), u;
    std::cout << "w= v + v = " << w << "\n"; 

    u= w - v;
    std::cout << "u= w - v = " << u << "\n"; 

    u= -v;
    std::cout << "u= -v = " << u << "\n"; 
    
    unroll<3>(u)= 4. * w - v;
    std::cout << "unroll<3>(u)= 4. * w - v = " << u << "\n"; 

    extra_test(u, v, u[0]);
}
 

int main(int, char**)
{
    using namespace mtl;
    using mtl::vec::parameters;

    dense_vector<float>   u(5);
    dense_vector<double>  x(5);
    dense_vector<std::complex<double> >  xc(5);

    std::cout << "Testing vector operations\n";

    test(u, "test float");
    test(x, "test double");

    test(xc, "test complex<double>");

    dense_vector<float, parameters<row_major> >   ur(5);
    test(ur, "test float in row vector");

    // dense_vector<dense_vector<float> > z(5, u), z2(5);
    // z2= 4.5f * z;
    // z2= z + z;

    return 0;
}
 














