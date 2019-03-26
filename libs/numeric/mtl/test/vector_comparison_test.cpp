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


#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  
    

template <typename VectorU, typename VectorV>
void test(VectorU u, VectorV v, const char* name)
{
    u= 4; v= 4;
    std::cout << "\n" << name << "\nu and v are " << (u == v ? "equal.\n" : "different.\n");
    MTL_THROW_IF(!(u == v), mtl::runtime_error("u and v should be equal."));
    
    std::cout << "u and v are " << (u != v ? "different.\n" : "equal.\n");
    MTL_THROW_IF(u != v, mtl::runtime_error("u and v should be equal."));

    // std::cout << "u is " << u << ", v is " << v << '\n';
    std::cout << "u and v are " << (u == 2*v ? "equal.\n" : "different.\n");
    MTL_THROW_IF(u == 2*v, mtl::runtime_error("u and v should be different."));

    std::cout << "u and v are " << (u != 2*v ? "different.\n" : "equal.\n");
    MTL_THROW_IF(!(u != 2*v), mtl::runtime_error("u and v should be different."));
}
 

int main(int, char**)
{
    using namespace mtl;
    using mtl::vec::parameters;

    dense_vector<int>     u(5);
    dense_vector<float>   v(5);
    dense_vector<double>  x(5);
    dense_vector<std::complex<double> >  xc(5);

    std::cout << "Testing vector comparison\n";

    test(u, u, "compare int with int");
    test(v, v, "compare float with float");
    test(x, x, "compare double with double");
    test(v, x, "compare float with double");

    test(xc, xc, "compare complex<double> with complex<double>");
    test(x, xc, "compare double with complex<double>");

    return 0;
}
