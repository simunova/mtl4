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
#include <string>

#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>

#include <boost/numeric/mtl/recursion/base_case_test.hpp>
#include <boost/numeric/mtl/recursion/base_case_matrix.hpp>
#include <boost/numeric/mtl/recursion/base_case_cast.hpp>


using namespace std;  
 
using mtl::recursion::base_case_matrix;
using mtl::recursion::simplify_base_case_matrix;


template <typename Matrix>
void test(Matrix& matrix)
{    
    using mtl::recursion::base_case_matrix;

    typedef mtl::recursion::max_dim_test_static<4>   base_test_type;
    base_test_type                              base_test;

    typedef typename base_case_matrix<Matrix, base_test_type>::type base_type;
    base_type base_matrix;
    cout << typeid(base_matrix).name() << "\n";
    
    Matrix sm= sub_matrix(matrix, 0, 4, 0, 4);
    // cout << typeid(simplify_base_case_matrix(sm, base_test)).name() << "\n";
    typename base_case_matrix<Matrix, base_test_type>::type simplified(simplify_base_case_matrix(sm, base_test));
}




int main(int, char**)
{
    using namespace mtl;

    typedef dense2D<int>                   d1t;
    typedef morton_dense<int, 0x55555553>  m1t; // col-major 4x4
    typedef morton_dense<int, 0x55555555>  m2t;
    typedef morton_dense<int, 0x5555555c>  m3t; // row-major 4x4
    typedef morton_dense<int, 0x555555f0>  m4t; // row-major 16x16

    d1t d1(8,8); m1t m1(8,8); m2t m2(8,8); m3t m3(8,8); m4t m4(8,8);

    test(d1);
    test(m1);
    test(m2);
    test(m3);
    test(m4);

    return 0;
} 
