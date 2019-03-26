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



using namespace mtl::recursion;
using namespace std;  


template <typename Matrix>
void test(Matrix& matrix, const char* name)
{
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    {
	mtl::mat::inserter<Matrix> ins(matrix);
	for (size_type i= 0; i < matrix.num_rows(); i++)
	    for (size_type j= 0; j < matrix.num_cols(); j++)
		if ((i + j) & 1)
		    ins(i, j) << i + 2*j;
    }

    // std::cout << "\n" << name << "\n" << matrix << "\n"; // Get us complaints from valgrind
    set_to_zero(matrix);
    std::cout << "\n" << name << "\n" << "should be empty now:\n" << matrix << "\n";
    typename Matrix::value_type zero(0.0);
    MTL_THROW_IF(matrix(0, 1) != zero, mtl::runtime_error("not properly set to zero"));
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

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(cdr, "Complex dense row major");
    test(md, "Morton N-order");
    test(d16r, "Hybrid 16 row-major");
    test(comp, "compressed2D");
    test(ccomp, "complex compressed2D");

    return 0;
}





