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
#include <boost/numeric/mtl/mtl.hpp>



using namespace std;  


template <typename MatrixA>
void test(MatrixA& A, unsigned dim1, unsigned dim2, const char* name)
{
    const unsigned max_print_size= 25;
    cout << "\n" << name << "\n";
    laplacian_setup(A, dim1, dim2);

    unsigned size= dim1 * dim2;
    mtl::dense_vector<double> v(size), b(size), r(size);

    for (unsigned i= 0; i < num_cols(A); i++)
	v[i]= A[12][i];
    b= 3.0;

    r= b + A * v;

    if (size <= max_print_size)
	cout << "A= \n" << A << "\n\nb + A * v = " << r << '\n';
    

    if (dim1 == 5 && dim2 == 5) {
	MTL_THROW_IF(abs(r[12] - 23.0) > 0.0001, mtl::runtime_error("r[12] should be 23.\n"));
    }

    r= b - A * v;

    if (size <= max_print_size)
	cout << "\n\nb - A * v = " << r << '\n';
    
    if (dim1 == 5 && dim2 == 5) {
	MTL_THROW_IF(abs(r[12] + 17.0) > 0.0001, mtl::runtime_error("r[12] should be -17.\n"));
    }

    // typedef mtl::dense_vector<double, mtl::parameters<mtl::row_major> >  vrt;
    // vrt rt, vt(trans(v));

    // rt= v * A;
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

    test(cr, dim1, dim2, "Row-major sparse");
    test(cc, dim1, dim2, "Column-major sparse");

    test(dr, dim1, dim2, "Row-major dense");
    test(dc, dim1, dim2, "Column-major dense");

    return 0;
}
 














