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
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  

template <typename Matrix, typename Indirect>
void check(const Indirect& I, const char* error)
{
    Matrix C(3, 2);
    C= 4, 2,
       2, 0,
       5, 3;

    C-= I;
    if (one_norm(C) > 0.01) throw error;
}

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    hessian_setup(A, 1.0);
    std::cout << "\n" << name << "\nA is: \n" << A;
    
    mtl::iset rows, cols;
    rows= 2, 0, 3;
    cols= 2, 0;

    cout << "rows = " << rows << ", cols = " << cols << "\n";
    cout << "A[rows][cols] is: \n" << A[rows][cols] << "\n";

    mtl::mat::indirect<Matrix> B(A[rows][cols]);
    cout << "B is\n" << B;
    check<Matrix>(B, "Wrong value after copy constructor");

    Matrix D(3, 2), E;
    D= B;
    check<Matrix>(D, "Wrong value after assignment");

    E= D + B;
    E/= 2;
    check<Matrix>(D, "Wrong value after addition");
}


int main(int, char**)
{
    using namespace mtl;
    const unsigned size= 5; 

    dense2D<double>                                  dc(size, size-2);
    dense2D<double, mat::parameters<col_major> >  dcc(size, size-2);
    dense2D<float>                                   fc(size, size-2);
    morton_dense<double,  morton_mask>               mdc(size, size-2);
    morton_dense<double, doppled_32_col_mask>        mcc(size, size-2);
    compressed2D<double>                             cc(size, size-2);
    compressed2D<double, mat::parameters<col_major> >  ccc(size, size-2);

    test(dc, "dense2D");
    test(dcc, "dense2D col-major");
    test(fc, "dense2D float");

    test(mdc, "pure Morton");
    test(mcc, "Hybrid col-major");
    test(cc, "Compressed");
    test(ccc, "Compresse col-majord");

    return 0;
}
