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

template <typename Vector>
void test(const char* name, const Vector& v)
{
    std::cout << name << ": v is " << v << ", trans(v) is " << trans(v) << ", trans(v) * v is " 
	      << trans(v) * v << "\n";
    MTL_THROW_IF(std::abs(std::abs(trans(v) * v) - 61.0) > 0.1, mtl::runtime_error("Wrong result"));
}

template <typename Vector>
void shape(const Vector&)
{
    std::cout << typeid(typename mtl::ashape::ashape<Vector>::type).name() << "\n";

}

template <typename Op1, typename Op2>
void shape(const Op1&, const Op2&)
{
    using namespace mtl;

    typedef typename ashape::mult_op<typename ashape::ashape<Op1>::type, typename ashape::ashape<Op2>::type>::type op_type;
    std::cout << "Op type is " << typeid(op_type).name() << "\n";

    typedef boost::is_same<op_type, ashape::rvec_cvec_mult>                same;
    typedef typename same::type                                            same_type;
    std::cout << "same type is " << typeid(same_type).name() << "\n";

    typedef boost::enable_if<same, float>          enable;
    typedef typename enable::type                  enable_type;
    std::cout << "enable type is " << typeid(enable_type).name() << "\n";

    typedef typename vec::detail::dot_result<Op1, Op2>::type res_type;
    std::cout << "res is " << typeid(res_type).name() << "\n";

}

template <typename Vector>
void test2(const char* name, const Vector& v)
{
    // shape(v); shape(trans(v));
    // shape(v, trans(v));
#if 1
    std::cout << name << ": v is " << v << ", trans(v) is " << trans(v) << ", v * trans(v) is " 
	      << v * trans(v) << "\n";
    MTL_THROW_IF(std::abs(std::abs(v * trans(v)) - 61.0) > 0.1, mtl::runtime_error("Wrong result"));
#endif 
}


int main(int, char**)
{

    typedef mtl::vec::fixed::dimension<3> fsize;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::row_major, fsize, true> >     rf; rf= 3., 4., 6.;
    mtl::dense_vector<float, mtl::vec::parameters<mtl::col_major, fsize, true> >     cf; cf= 3., 4., 6.;

    mtl::dense_vector<float, mtl::vec::parameters<mtl::row_major> >                  rd(3); rd= 3., 4., 6.;
    mtl::dense_vector<float>                                                            cd(3); cd= 3., 4., 6.;

    mtl::dense_vector<std::complex<double> >                                            rdc(3); rdc= 3., 4., 6.;

    mtl::dense2D<float>  A(3, 3);
    A= 2, 3, 4,
       3, 4, 6,
       7, 6, 9;

    test2("Row vector fixed size", rf);

    test("Column vector fixed size", cf);
    test2("Row vector", rd);
    test("Column vector", cd);
    //test2("Row vector complex", rdc);
#if 0
    test2("Matrix row", A[mtl::iall][1]);
    test("Matrix column", A[1][mtl::iall]);
#endif
    return 0;
}
