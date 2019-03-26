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
#include <cassert>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>

using namespace mtl; using std::cout;

template <typename Matrix, typename Memory>
void test(const char* text, const Matrix& A, Memory& m)
{
    // we can check some prerequisites
    assert(num_rows(A) <= num_rows(m));
    assert(num_cols(A) <= num_cols(m));
    assert(num_rows(m) == num_cols(m)); // m should be square


    // Define a tempory matrix with the size of A in the memory of m
    Memory tmp(num_rows(A), num_cols(A), &m[0][0]); 
    tmp= 0.0;

    // Define a bound (in our case of 16 x 16) and use it in the recursator
    int                          bound= num_rows(m);
    mat::recursator<Memory>    rec(tmp, bound);

    cout << "For matrix " << text << "\n";
    cout << "north_west(tmp)\n" << *north_west(rec);
    cout << "Is north_west empty or full: " << is_empty(north_west(rec)) 
	 << "  " << is_full(north_west(rec)) << "\n";
    cout << "south_west(tmp)\n" << *south_west(rec);
    cout << "Is south_west empty or full: " << is_empty(south_west(rec)) 
	 << "  " << is_full(south_west(rec)) << "\n";
}



int main(int argc, char* argv[])
{
    typedef dense2D<double>                                 matrix_type;
    //typedef morton_dense<double, recursion::morton_z_mask>  matrix_type;

    // Memory donator, can be used for all smaller matrices (not necessarily powers of 2)
    //    tested on dense2D and morton_dense
    matrix_type                                             mem(16, 16);         

    // Smaller matrices. Add more for testing
    matrix_type                                             A(8, 16), B(4, 8), C(4, 3);
    
    test("A", A, mem);
    test("B", B, mem);
    test("C", C, mem);

    return 0;
}

