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
    typedef typename mtl::Collection<Matrix>::value_type   value_type;
    
    value_type array[][3]= {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
    A= array;

    cout << "\n" << name << "\n" << "A =\n" << A;

    int indices[]= {1, 2, 0};
    mtl::mat::traits::permutation<>::type P= mtl::mat::permutation(indices);
    cout << "\nP =\n" << P;    

    Matrix A2( P * A );
    cout << "\nP * A =\n" << A2;

    MTL_THROW_IF(A2[1][2] != value_type(9.), mtl::runtime_error("Wrong value after row permutation!"));

    Matrix A3( A2 * trans(P) );
    cout << "\nA2 * trans(P) =\n" << A3;

    MTL_THROW_IF(A3[1][2] != value_type(7.), mtl::runtime_error("Wrong value after column permutation!"));
}


int main(int, char**)
{
    using namespace mtl;
    
    dense2D<double>                                      dr;
    dense2D<double, mat::parameters<col_major> >      dc;
    morton_dense<double, recursion::morton_z_mask>       mzd;
    morton_dense<double, recursion::doppled_2_row_mask>  d2r;
    compressed2D<double>                                 cr;
    compressed2D<double, mat::parameters<col_major> > cc;

    dense2D<complex<double> >                            drc;
    compressed2D<complex<double> >                       crc;

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");
    test(drc, "Dense row major complex");
    test(crc, "Compressed row major complex");

	
    return 0;
}
