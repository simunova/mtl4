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


template <typename T>
inline void check(const T& result, double exp)
{
    MTL_THROW_IF(std::abs(result - exp) > 0.0001, mtl::runtime_error("dot product wrong"));
}


template <typename VectorU, typename VectorV>
void test(VectorU& u, VectorV& v, const char* name)
{
    //using mtl::dot;
    typedef typename mtl::Collection<VectorU>::size_type  size_type;
    for (size_type i= 0; i < size(v); i++)
	u[i]= i+1, v[i]= i+1;

    
    mtl::io::tout << name << "\n dot(u, v) = " << dot(u, v) << "\n"; mtl::io::tout.flush();
    check(dot(u, v), 285.0);

    mtl::io::tout << " dot<2>(u, v) = " << mtl::dot<2>(u, v) << "\n"; mtl::io::tout.flush();
    check(mtl::dot<2>(u, v), 285.0);

    mtl::io::tout << " dot<6>(u, v) = " << mtl::dot<6>(u, v) << "\n"; mtl::io::tout.flush();
    check(mtl::dot<6>(u, v), 285.0);
}
 

int main(int ,char**)
{
    using mtl::vec::parameters;
    const int size= 9;

    mtl::dense_vector<float>   u(size), v(size), w(size);
    mtl::dense_vector<double>  x(size), y(size), z(size);
    mtl::dense_vector<std::complex<double> >  xc(size), yc(size), zc(size);

    mtl::io::tout << "Testing vector operations\n";

    test(u, v, "test float");
    test(x, y, "test double");
    test(u, x, "test float, double mixed");
    test(xc, yc, "test complex<double>");
    test(x, yc, "test complex<double>, double mixed");

    mtl::dense_vector<float, parameters<mtl::row_major> >   ur(size), vr(size), wr(size);
    test(ur, vr, "test float in row vector");
    
    // test(ur, v, wr, "test float in mixed vector (shouldn't work)"); 

    return 0;
}
 














