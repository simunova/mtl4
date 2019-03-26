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


template <typename Matrix>
void test(const char* A_string, const Matrix&)
{
    std::cout << A_string << " is " << mtl::static_num_rows<Matrix>::value << "x"
	      << mtl::static_num_cols<Matrix>::value << "\n";

    if (mtl::static_num_rows<Matrix>::value != 3)
	throw "Wrong number of rows";
    if (mtl::static_num_cols<Matrix>::value != 2)
	throw "Wrong number of columns";
    if (mtl::static_size<Matrix>::value != 6)
	throw "Wrong size";
}


int main(int, char**)
{
    using namespace mtl;

    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 2>, true> fmat_para;
    float ma[3][2]= {{2., 3.}, {4., 5.}, {6., 7.}};
    
    dense2D<float, fmat_para>        A_dense(ma);
    morton_dense<float, recursion::morton_z_mask, fmat_para>   A_morton(ma);
    compressed2D<float, fmat_para>   A_compressed(ma);

    test("Dense matrix", A_dense);
    test("Morton dense matrix", A_morton);
    test("Compressed sparse matrix", A_compressed);

    return 0;
}

