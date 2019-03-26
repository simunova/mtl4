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


#include <boost/numeric/mtl/mtl.hpp>


using namespace std;  

template <typename T>
T inline my_value(std::size_t i, T) { return T(i); }

std::complex<double> 
inline my_value(std::size_t i, std::complex<double>)
{ return std::complex<double>(i, i+1); }

template <typename VectorU>
void test(VectorU& u, const char* name)
{
    using std::abs;

    typedef typename mtl::Collection<VectorU>::value_type value_type;
    typedef typename mtl::Collection<VectorU>::size_type  size_type;
    for (size_type i= 0; i < size(u); i++)
	u[i]= my_value(i+1, value_type());

    value_type dot_cmp= dot(u, u);

    bool wrong= abs(unary_dot(u) - dot_cmp) > 0.001,
	 not_two_norm= abs(two_norm(u) * two_norm(u) - abs(unary_dot(u))) > 0.001;
    std::cout << name << ": u = " << u << "\n unary_dot(u) = " << unary_dot(u) << " is " 
	      << (wrong ? "not" : "") << "equal with dot(u, u)\n";     
    if (wrong) 
	throw "unary_dot product wrong";

    std::cout << " two_norm(u) = " << two_norm(u) << ", unary_dot(u) is " << (not_two_norm ? "not" : "") << "equal with square of it\n";
    if (not_two_norm)
	throw "different from square of two_norm";

    std::cout << " unary_dot<2>(u) = " << mtl::unary_dot<2>(u) << "\n"; std::cout.flush();
    if (abs(mtl::unary_dot<2>(u) - dot_cmp) > 0.001) 
	throw "unary_dot product wrong";

    std::cout << " unary_dot<6>(u) = " << mtl::unary_dot<6>(u) << "\n"; std::cout.flush();
    if (abs(mtl::unary_dot<6>(u) - dot_cmp) > 0.001) 
	throw "unary_dot product wrong";
}
 

int main(int ,char**)
{
    using mtl::vec::parameters;
    const int size= 9;

    mtl::dense_vector<float>   u(size);
    mtl::dense_vector<double>  x(size);
    mtl::dense_vector<std::complex<double> >  xc(size);

    std::cout << "Testing vector operations\n";

    test(u, "test float");
    test(x, "test double");
    test(xc, "test complex<double>");

    mtl::dense_vector<float, parameters<mtl::row_major> >   ur(size);
    test(ur, "test float in row vector");
    
    return 0;
}
 














