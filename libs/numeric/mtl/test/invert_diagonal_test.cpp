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


using namespace std;


template <typename Matrix>
void test(Matrix& A, const char* name)
{
    using mtl::Collection;
    A.change_dim(5, 5); A= 0.0;
    {
	mtl::mat::inserter<Matrix>   ins(A);
	ins[0][0] << 7; ins[1][1] << 8; ins[1][3] << 2; ins[1][4] << 3;
	ins[2][2] << 2; ins[3][3] << 4; ins[4][4] << 9;
     }
    
    cout << name << "\nA is set to \n" << A;
    invert_diagonal(A);
    
    cout << "\nAfter inverting the diagonal, it is: \n" << A;
    MTL_THROW_IF(std::abs(A[0][0] - 1.0 / 7.0) > 0.0001, mtl::runtime_error("Wrong value after inverting diagonal!"));
    MTL_THROW_IF(std::abs(A[1][3] - 2.0) > 0.0001, mtl::runtime_error("Wrong value after inverting diagonal!"));

}

int main(int, char**)
{
    mtl::dense2D<double>                                                dr;
    mtl::dense2D<double, mtl::mat::parameters<mtl::col_major> >      dc;
    mtl::morton_dense<double, mtl::recursion::morton_z_mask>            mzd;
    mtl::morton_dense<double, mtl::recursion::doppled_2_row_mask>       d2r;
    mtl::compressed2D<double>                                           cr;
    mtl::compressed2D<double, mtl::mat::parameters<mtl::col_major> > cc;

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");

    return 0;
}
