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

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/mtl/mtl.hpp> 


template <typename Vector, typename Ashape>
void test(const char* name, Vector& v, Ashape)
{
    v= 3, 4, 6;
    std::cout << name << ": v is " << v << ", typeid is " << typeid(v).name() << "\n";
    std::cout << "trans(v) is " << trans(v) << ", typeid is " << typeid(trans(v)).name() << "\n";
    BOOST_STATIC_ASSERT((boost::is_same<typename mtl::ashape::ashape<Vector>::type, Ashape>::value));
}

template <typename Vector, typename Ashape>
void test2(const char* name, const Vector& v, Ashape)
{
    std::cout << name << ": v is " << v << ", typeid is " << typeid(v).name() << "\n";
    std::cout << "trans(v) is " << trans(v) << ", typeid is " << typeid(trans(v)).name() << "\n";
    // trans(v)[0]= 2.3; // must not compile
    BOOST_STATIC_ASSERT((boost::is_same<typename mtl::ashape::ashape<Vector>::type, Ashape>::value));
}


int main(int, char**)
{
    using namespace mtl;

    typedef mtl::vec::fixed::dimension<3> fsize;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::row_major, fsize, true> >     rf;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::col_major, fsize, true> >     cf;

    mtl::dense_vector<float, mtl::vec::parameters<mtl::row_major> >                  rd(3);
    mtl::dense_vector<float>                                                         cd(3);

    mtl::dense_vector<std::complex<double> >                                         dc(3);

    mtl::dense2D<float>  A(3, 3);
    A= 2, 3, 4,
       3, 4, 6,
       7, 6, 9;

    test("Row vector fixed size", rf, ashape::rvec<ashape::scal>());
    test("Column vector fixed size", cf, ashape::cvec<ashape::scal>());
    test("Row vector", rd, ashape::rvec<ashape::scal>());
    test("Column vector", cd, ashape::cvec<ashape::scal>());
    test("Column vector complex", dc, ashape::cvec<ashape::scal>());
    test2("Matrix column", A[mtl::iall][1], ashape::cvec<ashape::scal>());
    test2("Matrix row", A[1][mtl::iall], ashape::rvec<ashape::scal>());

    return 0;
}
