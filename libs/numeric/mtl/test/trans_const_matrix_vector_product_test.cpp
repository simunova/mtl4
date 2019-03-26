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

// Test only compilability

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;


template <typename Matrix>
void test(Matrix& A, const char* name)
{
    A.change_dim(5, 5); A= 0.0;
    {
	mtl::mat::inserter<Matrix>   ins(A);
	ins[0][0] << 7; ins[1][1] << 8; ins[1][3] << 2; ins[1][4] << 3;
	ins[2][2] << 2; ins[3][3] << 4; ins[4][4] << 9;
    }
    
    double xa[] = {1, 2, 3, 4, 5};
    mtl::dense_vector<double> x(xa), b;
    
    const Matrix B(A);

    b= trans(A) * x;

    typedef mtl::transposed_view<Matrix>       trans_type;
    typedef mtl::transposed_view<const Matrix> ctrans_type;

#if 0
    cout << "Type of A" << typeid(typename mtl::ashape::ashape<Matrix>::type).name() << "\n";
    cout << "Type of B" << typeid(typename mtl::ashape::ashape<const Matrix>::type).name() << "\n";

    cout << "Type of trans(A)" << typeid(typename mtl::ashape::ashape<mtl::transposed_view<Matrix> >::type).name() << "\n";
    cout << "Type of trans(B)" << typeid(typename mtl::ashape::ashape<mtl::transposed_view<const Matrix> >::type).name() << "\n";
#endif

    cout << name <<":\n Type of enable_if_matrix<trans(A)> " << typeid(typename mtl::traits::enable_if_matrix<trans_type>::type).name() << "\n";
    cout << "Type of enable_if_matrix<trans(B)> " << typeid(typename mtl::traits::enable_if_matrix<ctrans_type>::type).name() << "\n";

    b= B * x;
    b= trans(B) * x;
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
#if 0
    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
#endif
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");

    return 0;
}
