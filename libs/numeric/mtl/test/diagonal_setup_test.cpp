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
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;  

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    typedef typename mtl::Collection<Matrix>::size_type     size_type;
    cout << "\n" << name << "\n";
    typename mtl::Collection<Matrix>::value_type four(4.0), one(1.0), zero(0.0);

    diagonal_setup(A, 1.0);
    cout << "Diagonal matrix:\n" << A << "\n";
    for (size_type r= 0; r < num_rows(A); ++r)
	for (size_type c= 0; c < num_cols(A); ++c)
	    if (r == c && A[r][c] != one) {throw "wrong diagonal";}
	    else MTL_THROW_IF(r != c && A[0][1] != zero, mtl::runtime_error("wrong off-diagonal"));

    A= 4.0;
    cout << "Diagonal matrix:\n" << A << "\n";
    for (size_type r= 0; r < num_rows(A); ++r)
	for (size_type c= 0; c < num_cols(A); ++c)
	    if (r == c && A[r][c] != four) {throw "wrong diagonal";}
	    else MTL_THROW_IF(r != c && A[0][1] != zero, mtl::runtime_error("wrong off-diagonal"));
}



int main(int, char**)
{
    using namespace mtl;
    unsigned size= 7; 

    dense2D<double>                                      dr(size, size), dr2(size, size+1), dr3(size+1, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);
    morton_dense<double, recursion::morton_z_mask>       mzd(size, size);
    morton_dense<double, recursion::doppled_2_row_mask>  d2r(size, size);
    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);

    dense2D<complex<double> >                            drc(size, size);
    compressed2D<complex<double> >                       crc(size, size);

    test(dr, "Dense row major");
    test(dr2, "Dense row major (one column more)");
    test(dr3, "Dense row major (one row more)");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");
    test(drc, "Dense row major complex");
    test(crc, "Compressed row major complex");

    return 0;
}
