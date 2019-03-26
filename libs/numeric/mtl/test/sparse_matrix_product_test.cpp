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
#include <boost/numeric/mtl/matrix/laplacian_setup.hpp> 
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>


using namespace std;  

template <typename MatrixA, typename MatrixB>
void test(MatrixA& a, MatrixB& b, unsigned dim1, unsigned dim2, const char* name)
{
    using mtl::Collection;

    const unsigned max_print_size= 25;
    cout << "\n" << name << "\n";
    laplacian_setup(a, dim1, dim2);
    laplacian_setup(b, dim1, dim2);

    unsigned size= dim1 * dim2;
    mtl::compressed2D<double>  c(size, size);
    c= a * b;

    if (size <= max_print_size)
	cout << "A= \n\n" << a << "B= \n\n" << b << "A*B= \n\n" << c << "\n";

    // Check for stencil below in the middle of the matrix
    //        1
    //     2 -8  2
    //  1 -8 20 -8  1
    //     2 -8  2
    //        1    
    if (dim1 == 5 && dim2 == 5) {
	typename Collection<MatrixA>::value_type twenty(20.0), two(2.0), one(1.0), 
	                                         zero(0.0), minus_eight(-8.0);
	MTL_THROW_IF(c[12][12] != twenty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(c[12][13] != minus_eight, mtl::runtime_error("wrong east neighbor"));
	MTL_THROW_IF(c[12][14] != one, mtl::runtime_error("wrong east east neighbor"));
	MTL_THROW_IF(c[12][15] != zero, mtl::runtime_error("wrong zero-element"));
	MTL_THROW_IF(c[12][17] != minus_eight, mtl::runtime_error("wrong south neighbor"));
	MTL_THROW_IF(c[12][18] != two, mtl::runtime_error("wrong south east neighbor"));
	MTL_THROW_IF(c[12][22] != one, mtl::runtime_error("wrong south south neighbor"));
    }

    c+= a * b;

    if (size <= max_print_size)
	cout << "C+= A*B= \n\n" << c << "\n";

    // Check for stencil, must be doubled now
    if (dim1 == 5 && dim2 == 5) {
	typename Collection<MatrixA>::value_type forty(40.0), four(4.0);
	MTL_THROW_IF(c[12][12] != forty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(c[12][18] != four, mtl::runtime_error("wrong south east neighbor"));
    }

    c-= a * b;

    if (size <= max_print_size)
	cout << "C-= A*B= \n\n" << c << "\n";

    // Check for stencil, must be A*B now
    if (dim1 == 5 && dim2 == 5) {
	typename Collection<MatrixA>::value_type twenty(20.0), two(2.0);
	MTL_THROW_IF(c[12][12] != twenty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(c[12][18] != two, mtl::runtime_error("wrong south east neighbor"));
    }
}



int main(int argc, char* argv[])
{
    using namespace mtl;

    unsigned dim1= 5, dim2= 5;

    if (argc > 2) {dim1= atoi(argv[1]);	dim2= atoi(argv[2]);}
    unsigned size= dim1 * dim2; 

    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);

    dense2D<complex<double> >                            drc(size, size);
    compressed2D<complex<double> >                       crc(size, size);

    test(cr, cr, dim1, dim2, "Row-major times row-major");
    test(cr, cc, dim1, dim2, "Row-major times column-major");

    test(cc, cr, dim1, dim2, "Column-major times row-major");
    test(cc, cc, dim1, dim2, "Column-major times column-major");

    return 0;
}
