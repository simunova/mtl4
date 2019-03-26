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

template <typename Matrix, typename Reorder>
void test_rows(const Matrix& A, const Reorder& r)
{
    typedef typename mtl::Collection<Matrix>::value_type   value_type;
    Matrix BB(reorder_matrix_rows(r, A));
    cout << "\nreorder_matrix_rows(A) =\n" << BB;
    MTL_THROW_IF(BB[1][0] != value_type(4.), mtl::runtime_error("Wrong value after row reordering!"));
}

template <typename Reorder, typename Value, typename Parameters>
void test_rows(const mtl::compressed2D<Value, Parameters>&, const Reorder&) {}


template <typename Matrix>
void test(Matrix& A, const char* name)
{
    typedef typename mtl::Collection<Matrix>::value_type   value_type;
    
    value_type array[][3]= {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}};
    A= array;

    cout << "\n" << name << "\n" << "A =\n" << A;

    int reordering[]= {2, 1};
    mtl::mat::traits::reorder<>::type  R= mtl::mat::reorder(reordering);
    cout << "\nR =\n" << R;    

    Matrix B(R * A);
    cout << "\nB= R * A =\n" << B;
    
    MTL_THROW_IF(B[1][0] != value_type(4.), mtl::runtime_error("Wrong value after row reordering!"));

    test_rows(A, reordering);
    Matrix B2(B * trans(R));
    cout << "\nB * trans(R) =\n" << B2;
    
    MTL_THROW_IF(B2[1][0] != value_type(6.), mtl::runtime_error("Wrong value after column reordering!"));    
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
