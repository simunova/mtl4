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

// #define MTL_VERBOSE_TEST

#include <iostream>
#include <algorithm>

#include <boost/numeric/mtl/mtl.hpp>


#ifdef MTL_WITH_MOVE

// Everything in the test is double
// Good enough for the moment

using namespace std;
using namespace mtl;
using mtl::io::tout;	


// Return a matrix with move semantics
// Return also the address of the first entry to be sure that it is really moved
template <typename Matrix>
Matrix f(const Matrix&, double*& a00)
{
    Matrix A(3, 3);
    A= 5.0;
    a00= &A.data[0];
    return A;
}

template <typename Matrix>
void print(const Matrix& matrix, const double* p)
{
    tout << "Data was " << (&matrix.data[0] == p ? "moved/shared.\n" : "copied.\n");
}

template <typename Matrix>
void test(const Matrix&, const char* text)
{
    tout << '\n' << text << '\n';

    double *p;
    Matrix A(3, 3);
    A= 0.0;
   
    tout << "A= f(A, p);\n";
    A= f(A, p);
    print(A, p);

    MTL_THROW_IF(A.data[0] != 5.0, mtl::runtime_error("Wrong value moving, should be 5.0!"));

    // Should be only on heap
    MTL_THROW_IF(!traits::is_static<Matrix>::value && &A.data[0] != p, mtl::runtime_error("Non-static matrix must be moved!"));
    MTL_THROW_IF(traits::is_static<Matrix>::value && &A.data[0] == p, mtl::runtime_error("Static matrix must be copied!"));

    tout << "Matrix B= f(A, p);\n";
    Matrix B(f(A, p));
    print(B, p);

    MTL_THROW_IF(B.data[0] != 5.0, mtl::runtime_error("Wrong value moving, should be 5.0!"));
    
    // Currently data are only moved on VS
    MTL_THROW_IF(!traits::is_static<Matrix>::value && &B.data[0] != p, mtl::runtime_error("Non-static matrix must be moved!"));
    // static data must be copied but that can be elided

    // This type is guarateed to be different to f's return type
    // In this case the matrix MUST be copied
    morton_dense<double, recursion::doppled_2_row_mask> C(3, 3);

    tout << "C= f(A, p);  // C and A have different types\n";
    C= f(A, p);
    print(C, p);

    MTL_THROW_IF(C.data[0] != 5.0, mtl::runtime_error("Wrong value trying to move, should be 5.0!"));
    MTL_THROW_IF(&C.data[0] == p, mtl::runtime_error("Matrix must be copied not moved!"));

    // Other matrix type, in this case the matrix MUST be copied
    morton_dense<double, recursion::morton_mask>   D(A);

    tout << "D(A);  // C and A have different types\n";
    print(D, &A.data[0]);

    MTL_THROW_IF(D.data[0] != 5.0, mtl::runtime_error("Wrong value in copy constructor, should be 5.0!"));
    MTL_THROW_IF(&D.data[0] == &A.data[0], mtl::runtime_error("Matrix must be copied not moved!"));

    p= &A.data[0];
    Matrix E(std::move(A));
    tout << "E(std::move(A));\n";
    print(E, p);
    // MTL_THROW_IF(&E.data[0] != p, mtl::runtime_error("Matrix must be moved not copied!"));
    MTL_THROW_IF(!traits::is_static<Matrix>::value && &E.data[0] != p, mtl::runtime_error("Non-static matrix must be moved!"));
    MTL_THROW_IF(traits::is_static<Matrix>::value && &E.data[0] == p, mtl::runtime_error("Static matrix must be copied!"));

}

template <typename Matrix>
void sub_matrix_test(const Matrix& A)
{
    Matrix F= sub_matrix(A, 0, 1, 0, 1);

    tout << "Matrix F= sub_matrix(A, 0, 1, 0, 1);\n";
    print(F, &A.data[0]);

    MTL_THROW_IF(&F.data[0] != &A.data[0], mtl::runtime_error("Sub-matrix must be referred to not copied!"));

    tout << "F= sub_matrix(A, 1, 2, 1, 2);\n";
    F= sub_matrix(A, 1, 2, 1, 2);    
    print(F, &A[1][1]);

    MTL_THROW_IF(&F.data[0] == &A[1][1], mtl::runtime_error("Matrix must be copied not referred to!"));

    Matrix G= clone(sub_matrix(A, 0, 1, 0, 1));

    tout << "Matrix G= clone(sub_matrix(A, 0, 1, 0, 1));\n";
    print(G, &A.data[0]);

    MTL_THROW_IF(&G.data[0] == &A.data[0], mtl::runtime_error("Sub-matrix must be forced to copy!"));
}



template <typename Matrix>
void dense_test(const Matrix& m, const char* text)
{
    test(m, text);
    sub_matrix_test(m);
}
#endif

int main(int, char*[])
{
#ifdef MTL_WITH_MOVE
    dense2D<double>                                   dr(3, 3);
    dense2D<double, mat::parameters<col_major> >      dc(3, 3);
    morton_dense<double, recursion::morton_z_mask>    mzd(3, 3);

    dense_test(dr, "Dense matrix");
    dense_test(dc, "Column-major dense matrix");
    dense_test(mzd, "Morton-order z-mask");

    // Check for data on stack
    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3>, true> fmat_para;
    dense2D<double, fmat_para>                        drs;
    test(drs, "Dense matrix on stack");

    compressed2D<double>                              crs(3, 3);
    compressed2D<double, mat::parameters<col_major> > ccs(3, 3);

    test(crs, "CRS");
    test(ccs, "CCS");
#else
    std::cout << "Test disabled due to lack of move semantics.";
#endif

    return 0;
}
