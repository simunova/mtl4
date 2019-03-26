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
#include <boost/numeric/mtl/vector/dense_vector.hpp> 
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>


using namespace std;  

template <typename MatrixA>
void test(MatrixA& A, unsigned dim1, unsigned dim2, const char* name)
{
    const unsigned max_print_size= 25;
    cout << "\n" << name << "\n";
    laplacian_setup(A, dim1, dim2);

    unsigned size= dim1 * dim2;
    mtl::dense_vector<double, mtl::vec::parameters<mtl::row_major> > v(size);
    unsigned r= size == 25 ? 12 : 3;
    for (unsigned i= 0; i < num_cols(A); i++)
	v[i]= A[r][i];

    // Resulting vector has same value type as matrix
    typedef typename mtl::Collection<MatrixA>::value_type rvalue_type;
    mtl::dense_vector<rvalue_type, mtl::vec::parameters<mtl::row_major> > w(size);

    w= v * A;

    if (size <= max_print_size)
	cout << "A= \n" << A << "\n\nv= " << v << "\n\nv * A= " << w << "\n";

    // Same test as in matrix product: resulting vector corresponds to column 12
    // Check for stencil below in the middle of the matrix
    //        1
    //     2 -8  2
    //  1 -8 20 -8  1
    //     2 -8  2
    //        1    
    if (dim1 == 5 && dim2 == 5) {
	rvalue_type twenty(20.0), two(2.0), one(1.0), zero(0.0), minus_eight(-8.0);
	MTL_THROW_IF(w[12] != twenty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(w[13] != minus_eight, mtl::runtime_error("wrong east neighbor"));
	MTL_THROW_IF(w[14] != one, mtl::runtime_error("wrong east east neighbor"));
	MTL_THROW_IF(w[15] != zero, mtl::runtime_error("wrong zero-element"));
	MTL_THROW_IF(w[17] != minus_eight, mtl::runtime_error("wrong south neighbor"));
	MTL_THROW_IF(w[18] != two, mtl::runtime_error("wrong south east neighbor"));
	MTL_THROW_IF(w[22] != one, mtl::runtime_error("wrong south south neighbor"));
    }

    w+= v * A;

    if (size <= max_print_size)
	cout << "w+= v*A= \n\n" << w << "\n";

    // Check for stencil, must be doubled now
    if (dim1 == 5 && dim2 == 5) {
	rvalue_type forty(40.0), four(4.0);
	MTL_THROW_IF(w[12] != forty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(w[18] != four, mtl::runtime_error("wrong south east neighbor"));
    }

    w-= v * A;

    if (size <= max_print_size)
	cout << "w-= v*A= \n\n" << w << "\n";

    // Check for stencil, must be v*A now
    if (dim1 == 5 && dim2 == 5) {
	rvalue_type twenty(20.0), two(2.0);
	MTL_THROW_IF(w[12] != twenty, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(w[18] != two, mtl::runtime_error("wrong south east neighbor"));
    }
}



int main(int, char**)
{
    using namespace mtl;

    unsigned dim1= 5, dim2= 5;
    unsigned size= dim1 * dim2; 

    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);

    dense2D<double>                                      dr(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);

    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<6, 6> > fmat_para;
    dense2D<double, fmat_para>                           drf;
   
    test(cr, dim1, dim2, "Row-major sparse");
    test(cc, dim1, dim2, "Column-major sparse");

    test(dr, dim1, dim2, "Row-major dense");
    test(dc, dim1, dim2, "Column-major dense");

    test(drf, 2, 3, "Row-major dense with static size 6x6");

    return 0;
}
