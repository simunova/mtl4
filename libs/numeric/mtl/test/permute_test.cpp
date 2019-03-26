// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/permute.hpp>

using namespace std;


int main(int, char**)
{
    using namespace mtl;

    dense_vector<int> v(4), p(4), ref(4);
    v = 7, 9, 10, 11;
    cout << "v = " << v << endl;
    p = 3, 1, 2, 0;
    cout << "p = " << p << endl;
    ref = 11, 9, 10, 7;
    
    dense_vector<int> vp = permute(p, v);
    cout << "vp = " << vp << endl;
    MTL_THROW_IF(vp != ref, mtl::runtime_error("Wrong value after permute"));
    
    dense_vector<int> v2 = reverse_permute(p, vp);
    cout << "v2 = " << v2 << endl;
    MTL_THROW_IF(v2 != v, mtl::runtime_error("Wrong value after reverse_permute"));
   
    return 0;
}
