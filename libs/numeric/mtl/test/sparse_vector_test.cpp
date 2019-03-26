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
#include <boost/numeric/mtl/vector/sparse_vector.hpp>


using namespace std;  


int main(int , char**)
{
#ifndef __PGI // trouble with zip_ref
    mtl::sparse_vector<float> v(100);
    std::cout << v << '\n';

    v.insert(3, 7.6f); v.insert(80, 3.9f); v.insert(40, 2.5f);
    std::cout << v << '\n';

    v[17]= 8.7f;
    v[19];
    v[18]+= 4.6f;
    std::cout << v << '\n';

    MTL_THROW_IF(!v.exists(18) || std::abs(v[18] - 4.6f) > 0.001f, logic_error("v[18] should exist and be 4.6"));
    MTL_THROW_IF(!v.exists(19) || v[19] != 0.0f, logic_error("v[19] should exist and be 0"));
    MTL_THROW_IF(v.exists(20), logic_error("v[20] should be empty"));
    
    v.sort_on_data();
    std::cout << v << '\n';
    // for
    MTL_THROW_IF(std::abs(v.entry(2).second) < std::abs(v.entry(3).second), logic_error("Vector not sorted on value magnitudes"));

    v.sort_on_indices();
    std::cout << v << '\n';
    MTL_THROW_IF(v.entry(2).first > v.entry(3).first, logic_error("Vector not sorted on indices"));

    std::cout << "two_norm(v) is " << two_norm(v) << '\n';
    std::cout << "infinity_norm(v) is " << infinity_norm(v) << '\n';
#endif

    return 0;
}
 














