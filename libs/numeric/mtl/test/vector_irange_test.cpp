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



const unsigned sz= 5;

inline float f(float x) { return x; }
inline double f(double x) { return x; }

inline std::complex<double> f(std::complex<double> x) 
{ 
    return std::complex<double>(real(x), real(x)+1.0); 
}


template <typename Vector>
void test(Vector& v, const char* name)
{
    using std::abs; using std::cout;
    using mtl::irange; using mtl::imax;

    typedef typename mtl::Collection<Vector>::value_type T;
    std::cout << "\n" << name << "\n";

    for (unsigned i= 0; i < size(v); ++i)
	v[i]= f(T(i));


    Vector w(v[irange(2, 4)]);
    w[1]= f(T(8));
    cout << "w == " << w << "\n";

    MTL_THROW_IF(w[0] != f(T(2)), mtl::runtime_error("Wrong value in w"));
    MTL_THROW_IF(size(w) != 2, mtl::runtime_error("Wrong size of w"));

    MTL_THROW_IF(v[3] != f(T(8)), mtl::runtime_error("Cannot change v via w (correctly)"));

    Vector u( v[irange(2, imax)] );
    cout << "u == " << u << "\n";
    cout << "v == " << v << "\n";
    irange r(2 , 4);

    MTL_THROW_IF(u[0] != f(T(2)), mtl::runtime_error("Wrong value in u"));
    MTL_THROW_IF(size(u) != sz-2, mtl::runtime_error("Wrong size of u"));

    cout << "v[irange(2, 4)] == " << v[r] << "\n";
    --r ;
    cout << "v[irange(2, 3)] == " << v[r] << "\n";

#if 0 // only for column vectors
    T beta= trans(v[irange(2, 4)]) * w;
    cout << "trans(v) * w = " << beta << "\n";
#endif
}



int main(int, char**)
{
    using namespace mtl;
    dense_vector<float>                                                 cf(sz, 1.0);
    dense_vector<double>                                                cd(sz, 1.0);
    dense_vector<std::complex<double> >                                 cc(sz, 1.0);
    dense_vector<float, mtl::vec::parameters<row_major> >               rf(sz, 1.0);

    test(cf, "dense_vector<float>");
    test(cd, "dense_vector<double>");
    test(cc, "dense_vector<std::complex<double> >");
    test(rf, "dense_vector<float, parameters<row_major> >");

    return 0;
}
