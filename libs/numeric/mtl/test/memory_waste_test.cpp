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

int main(int, char**)
{
    using namespace std;
    
    const unsigned s= 1000;
    mtl::compressed2D<double> A(s, s);
    A= 1.0;
    
    MTL_THROW_IF(A.data.capacity() > 2 * A.data.size(), mtl::runtime_error("Capacity too large"));
    MTL_THROW_IF(A.ref_major().capacity() > 2 * A.ref_major().size(), mtl::runtime_error("Capacity too large"));
    MTL_THROW_IF(A.ref_minor().capacity() > 2 * A.ref_minor().size(), mtl::runtime_error("Capacity too large"));

    return 0;
}
