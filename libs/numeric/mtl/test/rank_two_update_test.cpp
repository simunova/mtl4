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
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/laplacian_setup.hpp> 
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/rank_two_update.hpp>


using namespace std;  

inline double value(double)
{
    return 1.0;
}

inline complex<double> value(complex<double>)
{
    return complex<double>(1.0, 1.0);
}

inline double test_value(double)
{
    return 10.0;
}

inline complex<double> test_value(complex<double>)
{
    return complex<double>(10.0, -10.0);
}


template <typename Matrix>
void test(Matrix& matrix, const char* name)
{
    using mtl::conj; using mtl::Collection;
    const unsigned max_print_size= 25;

    cout << "\n" << name << "\n";
    set_to_zero(matrix);

    typename Collection<Matrix>::size_type          nr= num_rows(matrix), nc= num_cols(matrix);
    typedef typename Collection<Matrix>::value_type value_type;
    value_type                                      zero(0.0);
    mtl::dense_vector<value_type>                   x(nr, zero), y(nc, zero);

    x[1]= 1.0; x[2]= 2.0;
    value_type ref(0), v= value(ref);

    y[4]= 4.0*v; y[5]= 5.0*v; y[6]= 6.0*v;

    rank_two_update(matrix, x, y);
    if (nr <= max_print_size)
	cout << "\nx= " << x << "y= " << y << "matrix = \n" << matrix << "\n";

    MTL_THROW_IF(matrix[2][5] != test_value(v), mtl::runtime_error("wrong value"));
    MTL_THROW_IF(matrix[5][2] != conj(test_value(v)), mtl::runtime_error("wrong value"));
}



int main(int argc, char* argv[])
{
    using namespace mtl;

    cout << "matrix size must be at least 7 x 7\n";
    unsigned size= 7;
    if (argc > 1) size= atoi(argv[1]);     

    dense2D<double>                                      dr(size, size);
    dense2D<double, mat::parameters<col_major> >      dc(size, size);
    morton_dense<double, recursion::morton_z_mask>       mzd(size, size);
    morton_dense<double, recursion::doppled_2_row_mask>  d2r(size, size);
    compressed2D<double>                                 cr(size, size);
    compressed2D<double, mat::parameters<col_major> > cc(size, size);

    dense2D<complex<double> >                            drc(size, size);
    compressed2D<complex<double> >                       crc(size, size);

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
