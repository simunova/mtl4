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


template <typename Vector>
void test(const char* name, const Vector&)
{
    std::cout << name << " is " << mtl::static_num_rows<Vector>::value << "x"
	      << mtl::static_num_cols<Vector>::value << "\n";

    if (mtl::static_num_rows<Vector>::value != (mtl::traits::is_row_major<Vector>::value ? 1 : 2))
	throw "Wrong number of rows";
    if (mtl::static_num_cols<Vector>::value != (mtl::traits::is_row_major<Vector>::value ? 2 : 1))
	throw "Wrong number of columns";
    if (mtl::static_size<Vector>::value != 2)
	throw "Wrong size";
}


int main(int, char**)
{
    using namespace mtl;

    typedef mtl::vec::parameters<tag::col_major, mtl::vec::fixed::dimension<2>, true> col_para;
    typedef mtl::vec::parameters<tag::row_major, mtl::vec::fixed::dimension<2>, true> row_para;
    float va[2]= {3., 4.};
    dense_vector<float, col_para>   v_col(va);
    dense_vector<float, row_para>   v_row(va);

    test("Dense column vector", v_col);
    test("Dense row vector", v_row);

    return 0;
}

