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
#include <complex>
#include <cmath>

#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/transposed_view.hpp>

#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/copy.hpp>



using namespace mtl::recursion;
using namespace std;  


template <typename Matrix>
void init_matrix(Matrix& matrix, int offset= 0)
{
    set_to_zero(matrix);
    mtl::mat::inserter<Matrix> ins(matrix);
    for (unsigned i= 0; i < matrix.num_rows(); i++)
	for (unsigned j= 0; j < matrix.num_cols(); j++)
	    if ((i + j + offset) & 1)
		ins(i, j) << i + 2*j;
}


template <typename MatrixSrc, typename MatrixDest>
void test(MatrixSrc& src, const char* name_src, 
	  MatrixDest&, const char* name_dest)
{
    std::cout << "\nAssign " << name_src << " to " << name_dest << "\n";

    MatrixDest new_dest(3, 4);
    init_matrix(new_dest);
    std::cout << "dest initialized:\n" << new_dest << "\n";

    init_matrix(src, 1);
    std::cout << "source initialized:\n" << src << "\n";

    // mat::copy(src, new_dest);
    new_dest.change_dim(num_rows(src), num_cols(src));
    new_dest= src;
    std::cout << "dest after assignment:\n" << new_dest << "\n\n";

    MTL_THROW_IF(new_dest.num_rows() != 5 || new_dest.num_cols() != 7, mtl::runtime_error("wrong dimension"));
    typename MatrixDest::value_type zero(0.0), three(3.0);
    MTL_THROW_IF(new_dest(1, 2) != zero, mtl::runtime_error("m[1][2] should be zero"));
    MTL_THROW_IF(new_dest(1, 1) != three, mtl::runtime_error("m[1][1] should be three"));

}


int main(int, char**)
{
    using namespace mtl;

    dense2D<double>                                dr(5, 7);
    dense2D<double, mat::parameters<col_major> > dc(5, 7);
    dense2D<std::complex<double> >                 cdr(5, 7);
    morton_dense<double,  morton_mask>             md(5, 7);
    morton_dense<double,  doppled_16_row_mask>     d16r(5, 7);
    compressed2D<double>                           comp(5, 7);
    compressed2D<std::complex<double> >            ccomp(5, 7);

    test(dc, "Dense column major", dc, "Dense column major");
    test(dc, "Dense column major", cdr, "Complex dense row major");
    test(dc, "Dense column major", md, "Morton N-order");
    test(dc, "Dense column major", d16r, "Hybrid 16 row-major");
    test(dc, "Dense column major", comp, "compressed2D");
    test(dc, "Dense column major", ccomp, "complex compressed2D");

    test(d16r, "Hybrid 16 row-major", dr, "Dense row major");
    test(d16r, "Hybrid 16 row-major", dc, "Dense column major");
    test(d16r, "Hybrid 16 row-major", cdr, "Complex dense row major");
    test(d16r, "Hybrid 16 row-major", md, "Morton N-order");
    test(d16r, "Hybrid 16 row-major", d16r, "Hybrid 16 row-major");
    test(d16r, "Hybrid 16 row-major", comp, "compressed2D");
    test(d16r, "Hybrid 16 row-major", ccomp, "complex compressed2D");

    test(comp, "compressed2D", dr, "Dense row major");
    test(comp, "compressed2D", dc, "Dense column major");
    test(comp, "compressed2D", cdr, "Complex dense row major");
    test(comp, "compressed2D", md, "Morton N-order");
    test(comp, "compressed2D", d16r, "Hybrid 16 row-major");
    test(comp, "compressed2D", comp, "compressed2D");
    test(comp, "compressed2D", ccomp, "complex compressed2D");

    return 0;
}


