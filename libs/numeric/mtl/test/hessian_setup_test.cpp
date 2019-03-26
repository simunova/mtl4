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
void test(MatrixA& A, const char* name)
{
    hessian_setup(A, 1.0);
    std::cout << "Hessian setup for " << name << "\n";
    if (num_rows(A) < 10)
	std::cout << "A is\n" << A;

    MTL_THROW_IF(num_rows(A) >= 4 && std::abs(A[3][2] - 5) > 0.001, mtl::runtime_error("A[3][2] should be 5"));
}
 

int main(int argc, char* argv[])
{
    using namespace mtl;

    unsigned size= 5; 
    if (argc > 1) size= atoi(argv[1]); 
    if (size < 2) size= 2;

    dense2D<double>                                       da(size, size-1); 
    dense2D<double, mat::parameters<col_major> >       dca(size, size-1);
    dense2D<float>                                        fa(size, size-1);

    morton_dense<double,  morton_mask>                    mda(size, size-1);
    morton_dense<double, doppled_32_col_mask>             mca(size, size-1);
    morton_dense<double, doppled_32_row_mask>             mra(size, size-1);

    compressed2D<double>                                  cda(size, size-1); 
    compressed2D<double, mat::parameters<col_major> >  cdca(size, size-1);
    

    test(da, "dense2D");
    test(dca, "dense2D col-major");
    test(fa, "dense2D float");

    test(mda, "pure Morton");
    test(mca, "Hybrid col-major");
    test(mra, "Hybrid row-major");

    test(cda, "compressed2D");
    test(cdca, "compressed2D col-major");

    return 0;
}
 














