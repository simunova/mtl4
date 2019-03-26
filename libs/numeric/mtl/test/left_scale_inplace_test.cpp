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
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/dense2D.hpp> 
#include <boost/numeric/mtl/matrix/laplacian_setup.hpp> 
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/operation/left_scale_inplace.hpp>


using namespace std;  

template <typename MatrixA, typename MatrixB>
void test(MatrixA& a, MatrixB& b, unsigned dim1, unsigned dim2, const char* name)
{
    unsigned size= dim1 * dim2;
    MTL_THROW_IF(size == 0, mtl::runtime_error("Matrix size must be larger than 0 to make the test meaningful."));

    const unsigned max_print_size= 25;
    cout << "\n" << name << "\n";
    laplacian_setup(a, dim1, dim2);
    laplacian_setup(b, dim1, dim2);

    left_scale_inplace(2.0, a);
    if (size <= max_print_size)
	cout << "A= \n\n" << a << "\n";

    typename mtl::Collection<MatrixA>::value_type eight(8.0);
    MTL_THROW_IF(a[0][0] != eight, mtl::runtime_error("Scaling with scalar wrong"));

    left_scale_inplace(0.5, a);
    left_scale_inplace(b, a);

    if (size <= max_print_size)
	cout << "A= \n\n" << a << "B= \n\n" << b << "\n";

    // Check for stencil below in the middle of the matrix
    //        1
    //     2 -8  2
    //  1 -8 20 -8  1
    //     2 -8  2
    //        1    
    if (dim1 == 5 && dim2 == 5) {
	typename mtl::Collection<MatrixA>::value_type twenty(20.0), two(2.0), one(1.0), 
	                                              zero(0.0), minus_eight(-8.0);
	MTL_THROW_IF(a[12][12] != twenty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(a[12][13] != minus_eight, mtl::runtime_error("wrong east neighbor"));
	MTL_THROW_IF(a[12][14] != one, mtl::runtime_error("wrong east east neighbor"));
	MTL_THROW_IF(a[12][15] != zero, mtl::runtime_error("wrong zero-element"));
	MTL_THROW_IF(a[12][17] != minus_eight, mtl::runtime_error("wrong south neighbor"));
	MTL_THROW_IF(a[12][18] != two, mtl::runtime_error("wrong south east neighbor"));
	MTL_THROW_IF(a[12][22] != one, mtl::runtime_error("wrong south south neighbor"));
    }
}



int main(int argc, char* argv[])
{
    using namespace mtl;
    unsigned dim1= 5, dim2= 5;

    if (argc > 2) {dim1= atoi(argv[1]);dim2= atoi(argv[2]);}
    unsigned size= dim1 * dim2; 

    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);

    dense2D<double>                                      dr(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);

    test(cr, dr, dim1, dim2, "Row-major sparse scaled with row-major dense");
    test(cr, dc, dim1, dim2, "Row-major sparse scaled with column-major dense");
    test(cc, dr, dim1, dim2, "Column-major sparse scaled with row-major dense");
    test(cc, dc, dim1, dim2, "Column-major sparse scaled with column-major dense");

    test(dr, cr, dim1, dim2, "Row-major dense scaled with row-major sparse");
    test(dr, cc, dim1, dim2, "Row-major dense scaled with column-major sparse");
    test(dc, cr, dim1, dim2, "Column-major dense scaled with row-major sparse");
    test(dc, cc, dim1, dim2, "Column-major dense scaled with column-major sparse");

    return 0;
}
