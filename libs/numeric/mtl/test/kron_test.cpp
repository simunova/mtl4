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

int main()
{
    using namespace mtl;
    
    dense2D<int> A(2, 2), B(2, 2), C(4, 4);
    
    for (size_t r= 0; r < 2; ++r)
        for (size_t c= 0; c < 2; ++c) {
            A[r][c]= (r+1) * 10 + c+1;
            B[r][c]= (r+1) * 1000 + (c+1) * 100;
        }
        
    C= kron(A, B);
    std::cout << "kron(A, B) is\n" << C;
    
    MTL_THROW_IF(C[0][0] != 12100, mtl::runtime_error("Wrong value in C[0][0]"));
    MTL_THROW_IF(C[3][3] != 48400, mtl::runtime_error("Wrong value in C[3][3]"));

    return 0;
}
