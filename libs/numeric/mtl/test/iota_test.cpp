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
#include <vector>

#include <boost/numeric/mtl/mtl.hpp>

const unsigned sz= 5;

inline float f(float x) { return x; }
inline double f(double x) { return x; }

template <typename Vector>
void test(Vector& v, const char* name)
{ 
    typedef typename mtl::Collection<Vector>::value_type T;
    std::cout << "\n" << name << "\n";

    Vector w(v);
    iota(w);
    MTL_THROW_IF(w[4] != f(T(4)), mtl::runtime_error("Wrong value in w"));
    std::cout << "w == " << w << "\n";
    
    iota(w, T(5));
    MTL_THROW_IF(w[4] != f(T(9)), mtl::runtime_error("Wrong value in w"));
    std::cout << "w == " << w << "\n";

}



int main(int, char**)
{
    using namespace mtl;
    dense_vector<float>                                                 cf(sz, 0.0);
    dense_vector<double>                                                cd(sz, 0.0);
    dense_vector<float, mtl::vec::parameters<row_major> >            rf(sz, 0.0);

    test(cf, "dense_vector<float>");
    test(cd, "dense_vector<double>");
    test(rf, "dense_vector<float, parameters<row_major> >");

    return 0;
}
