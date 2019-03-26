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
#include <string>
#include <vector>
#include <complex>
#include <boost/timer.hpp>


// Define only if you really have an Opteron!!!!!!!!!!!!!!!
// #define MTL_USE_OPTERON_OPTIMIZATION

#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;
using namespace mtl::recursion; 
using namespace std;  
using assign::plus_sum; using assign::assign_sum; 

typedef bound_test_static<32>  test32_t;

// Use tiles of 2x4, usually works best on odin
// Any other tile sizes NxM can be tested with N * M <= 16 
typedef gen_tiling_dmat_dmat_mult_t<2, 4, plus_sum>    tiling_mult_t;

// Recursive multiplication with 32x32 blocks and tiling in the blocks
typedef gen_recursive_dmat_dmat_mult_t<tiling_mult_t, test32_t>  recursive_mult_t;

// Assembler code if applicable, only implemented on certain bit masks 
//   and only for double values (otherwise falls back on canonical impl. (slow))
typedef gen_platform_dmat_dmat_mult_t<plus_sum>     platform_mult_t;

// Recursive assembler version
typedef gen_recursive_dmat_dmat_mult_t<platform_mult_t, test32_t>  recursive_asm_mult_t;


int main(int argc, char* argv[])
{
    const int                                      n= 10;
    morton_dense<double,  doppled_32_row_mask>     A1(n, n), C1(n, n);
    morton_dense<double,  doppled_32_col_mask>     B1(n, n);
    recursive_mult_t                               m1;

    // fill matrix, e.g.
    hessian_setup(A1, 1.0);
    hessian_setup(B1, 1.0);

    // Recursive with fast 32x32 base case in C++
    m1(A1, B1, C1);
    cout << "Result with C++ is:\n" << C1 << "\n";

    

    // Now assembly code with shark teeth

    morton_dense<double,  shark_32_row_mask>     A2(n, n), C2(n, n);
    morton_dense<double,  shark_32_col_mask>     B2(n, n);
    recursive_asm_mult_t                         m2;

    // fill matrix, e.g.
    hessian_setup(A2, 1.0);
    hessian_setup(B2, 1.0);

    // Recursive with fast 32x32 base case in assembly
    m2(A2, B2, C2);
    cout << "Result with assembly is:\n" << C2 << "\n";

    return 0; 
}
