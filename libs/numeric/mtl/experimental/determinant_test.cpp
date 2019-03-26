// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/timer.hpp>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;

template <typename Matrix>
typename mtl::Collection<Matrix>::value_type
inline determinant(const Matrix& A, double eps= 0)
{
    mtl::dense_vector<std::size_t> P(num_rows(A)); 
    mtl::dense2D<typename mtl::Collection<Matrix>::value_type>   LU(A);  // A is copied in a dense matrix because lu overwrites its argument
    lu(LU, P, eps);
    return product(diagonal(LU));
}

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    laplacian_setup(A, 3, 4);
    // boost::timer t;
    cout << "\n" << name << ", determinant is " << determinant(A) << "\n";
    // cout << "It took " << t.elapsed() << "s\n";
}


int main(int, char**) 
{
    using namespace mtl;
    dense2D<double>                                      dr;
    dense2D<complex<double> >                            dz;
    dense2D<double, mat::parameters<col_major> >      dc;
    compressed2D<double>                                 cr;

    test(dr, "Row-major dense");
    test(dz, "Row-major dense with complex numbers");
    test(dc, "Column-major dense");
    test(cr, "Row-major compressed");
    
    return 0;
}
