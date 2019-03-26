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
#include <boost/numeric/mtl/mtl.hpp>

int main(int , char**)
{
    using namespace std;
    
    mtl::iset indices;

    std::cout << "indices = " << indices << '\n';
    MTL_THROW_IF(indices.size() != 0, mtl::runtime_error("Set should be empty."));

    indices.push_back(3); indices.push_back(5); 
    std::cout << "indices = " << indices << '\n';
    MTL_THROW_IF(indices.size() != 2, mtl::runtime_error("Set should have size 2."));
    
    mtl::iset i2(indices);
    std::cout << "Index sets are" << (indices == i2 ? "" : " not") << " equal\n";
    std::cout << "Index sets are" << (indices != i2 ? " not" : "") << " equal\n";

    std::size_t index_array[]= {3, 5, 2, 0};
    mtl::iset i3(index_array);
    std::cout << "i3 = " << i3 << '\n';

    i2= 3, 5, 2, 0;
    std::cout << "i2 = " << i2 << '\n';

    return 0;
}
