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
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/compressed2D.hpp> 
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/hermitian_view.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/mtl/operation/scale.hpp>
#include <boost/numeric/mtl/operation/hermitian.hpp>


using std::cout; using std::complex;  

typedef complex<double> ct;

double value(double)
{
    return 7.0;
}

complex<double> value(complex<double>)
{
    return ct(7.0, 1.0);
}



template <typename MatrixA, typename MatrixB, typename MatrixC>
void test(MatrixA&, MatrixB&, MatrixC&, const char* name)
{
    MatrixA a(7, 7); 
    MatrixB b(7, 7); 
    MatrixC c(7, 7);

    set_to_zero(a); 
    {
	typename MatrixA::value_type ref(0);
	mtl::mat::inserter<MatrixA>  ins(a);
	ins(2, 3) << value(ref);
	ins(4, 3) << value(ref) + 1.0;
	ins(2, 5) << value(ref) + 2.0;
    }
    std::complex<double> sum= a[4][3], diff= a[4][3];


    set_to_zero(b);
    {
	typename MatrixB::value_type ref(0);
	mtl::mat::inserter<MatrixB>  ins(b);
	ins(2, 2) << value(ref) + 3.0;
	ins(4, 3) << value(ref) + 4.0;
    }
    sum+= b[4][3];
    diff-= b[4][3];

    cout << "\n\n" << name << "\n";
    cout << "Original matrices:\nA=\n" << a << "B=\n" << b << "\n";

    c= a + b;
    
    cout << "C= A + B\n" << c << "\n";
    MTL_THROW_IF(c[4][3] != sum, mtl::runtime_error("wrong sum"));

    c= a + b + a + b;
    
    cout << "C= A + B + A + B\n" << c << "\n";
    MTL_THROW_IF(c[4][3] != 2.0*sum, mtl::runtime_error("wrong sum"));
    c+= a + b;
    
    cout << "C+= A + B\n" << c << "\n";
    MTL_THROW_IF(c[4][3] != 3.0*sum, mtl::runtime_error("wrong increment by sum"));

    cout << "A + B\n" << a+b << "\n";
    cout << "(A + B)[4][3] = " << (a+b)[4][3] << "\n\n";

    c-= a + b;
    cout << "C-= A + B\n" << c << "\n";
    MTL_THROW_IF(c[4][3] != 2.0*sum, mtl::runtime_error("wrong decrement by sum"));

    c= a - b;
    
    cout << "C= A - B\n" << c << "\n";
    MTL_THROW_IF(c[4][3] != diff, mtl::runtime_error("wrong difference"));

    c-= a - b;    
    cout << "C-= A - B\n" << c << "\n";
    MTL_THROW_IF(c[4][3] != 0.0, mtl::runtime_error("wrong decrement by difference"));
}



int main(int argc, char* argv[])
{
    using namespace mtl;
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


    test(dr, dr, dr, "Dense row major");
    test(dc, dr, dr, "Dense column major as sum of dense rows");
    test(dc, dr, dc, "Dense column major as sum of dense rows and column");

    test(mzd, mzd, mzd, "Morton Z-order");
    test(d2r, mzd, d2r, "Hybrid 2 row-major + Morton Z-order");

    test(cr, cr, cr, "Compressed row major");
    test(cc, cr, cc, "Compressed column major + row");

    test(drc, drc, drc, "Dense row major complex");
    test(drc, dc, drc, "Dense row major complex + column double");
    test(crc, crc, crc, "Compressed row major complex");
    test(crc, dc, crc, "Compressed row major complex + dense column major double");

    return 0;
}
