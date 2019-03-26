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
#include <vector>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;
// 

template <typename Vector>
void test(Vector& v, const char* name)
{
    using mtl::size; using mtl::num_rows; using mtl::num_cols;

    cout << "\n" << name << "\n";

    cout << "size(v) = " << size(v) << "\n";

    // Value is less critical main purpose of the test is to check compilibilit
    MTL_THROW_IF(size(v) != 3, mtl::runtime_error("Vector size should be 3"));

    cout << "num_rows(v) = " << num_rows(v) << "\n";

    // Value is less critical main purpose of the test is to check compilibilit
    MTL_THROW_IF(num_rows(v) != 3, mtl::runtime_error("Vector number of rows should be 3"));

    cout << "num_cols(v) = " << num_cols(v) << "\n";

    // Value is less critical main purpose of the test is to check compilibilit
    MTL_THROW_IF(num_cols(v) != 1, mtl::runtime_error("Vector number of columns should be 1"));

}


int main(int, char**)
{
    mtl::dense_vector<int>   mv(3, 3);
    std::vector<int>         sv(3, 3);
    int                      array[3];

    test(mv, "MTL dense_vector");
    test(sv, "STL vector");
    test(array, "int[3]");
	
    return 0;
}
