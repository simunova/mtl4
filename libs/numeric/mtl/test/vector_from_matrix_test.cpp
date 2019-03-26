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
    
template <typename Vector, typename Value>
void check(const char* name, const Vector& v, bool row_major, unsigned s, Value c0, Value c1)
{
    std::cout << name << " = " << v << "\n";

    MTL_THROW_IF(mtl::traits::is_row_major<Vector>::value != row_major, mtl::runtime_error("wrong orientation"));
    MTL_THROW_IF(size(v) != s, mtl::runtime_error("wrong size"));
    MTL_THROW_IF(v[0] != c0, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(v[1] != c1, mtl::runtime_error("wrong value"));
}



template <typename Matrix>
void test(Matrix& A, const char* name)
{
    using mtl::irange; using mtl::iall; using mtl::imax;
    typedef typename mtl::Collection<Matrix>::value_type value_type;
    typedef mtl::dense_vector<float>                     Vector;

    std::cout << "\n" << name << "  --- A = \n" << A; 

    Vector vc(3), vr(3);
    vc= 2, 5, 8; vr= 1, 2, 3;

    cout << "A[iall][1] == " << A[iall][1] << "\n";
    MTL_THROW_IF(one_norm(Vector(vc - A[iall][1])) > 0.1, mtl::runtime_error("Wrong column vector"));
    check("A[iall][1]", A[iall][1], false, 3, value_type(2), value_type(5));

    cout << "A[irange(1, imax)][1] == " << A[irange(1, imax)][1] << "\n";
    MTL_THROW_IF(one_norm(Vector(vc[irange(1, imax)] - A[irange(1, imax)][1])) > 0.1, mtl::runtime_error("Wrong column cub-vector"));
    check("A[irange(1,3)][1]", A[irange(1,3)][1], false, 2, value_type(5), value_type(8));

    typename mtl::ColumnInMatrix<Matrix>::type c(A[irange(1, imax)][1]);
    c[1]= 8.5;
    MTL_THROW_IF(A[2][1] != 8.5, mtl::runtime_error("Matrix modification (in column) did not work"));
    check("c= A[irange(1, imax)][1]", c, false, 2, value_type(5), value_type(8.5));

    check("A[1][iall]", A[1][iall], true, 3, value_type(4), value_type(5));
    check("A[1][irange(1,3)]", A[1][irange(1,3)], true, 2, value_type(5), value_type(6));

    typename mtl::RowInMatrix<Matrix>::type r(A[1][irange(1, imax)]);
    r[1]= 6.5;
    MTL_THROW_IF(A[1][2] != 6.5, mtl::runtime_error("Matrix modification (in row) did not work"));
    check("r= A[1][irange(1, imax)]", r, true, 2, value_type(5), value_type(6.5));
}
 
template <typename Matrix>
void test(const Matrix& A, const char* name)
{
    using mtl::irange; using mtl::iall; using mtl::imax;
    typedef typename mtl::Collection<Matrix>::value_type value_type;

    std::cout << "\n" << name << "  --- A = \n" << A; 
    check("A[iall][1]", A[iall][1], false, 3, value_type(2), value_type(5));
    check("A[irange(1,3)][1]", A[irange(1,3)][1], false, 2, value_type(5), value_type(8.5));

    typename mtl::ColumnInMatrix<const Matrix>::type c(A[irange(1, imax)][1]);
    check("c= A[irange(1, imax)][1]", c, false, 2, value_type(5), value_type(8.5));

    check("A[1][iall]", A[1][iall], true, 3, value_type(4), value_type(5));
    check("A[1][irange(1,3)]", A[1][irange(1,3)], true, 2, value_type(5), value_type(6.5));

    typename mtl::RowInMatrix<const Matrix>::type r(A[1][irange(1, imax)]);
    check("r= A[1][irange(1, imax)]", r, true, 2, value_type(5), value_type(6.5));
}


int main(int, char**)
{
    using namespace mtl;

    dense2D<float> A(3, 3);
    A= 1, 2, 3,
       4, 5, 6, 
       7, 8, 9;

    dense2D<float, mat::parameters<col_major> > B(A);
    
    test(A, "Row-major matrix");     // A changed !!!
    test(B, "Column-major matrix");  // B changed !!!

    const dense2D<float>                                 C(A);
    const dense2D<float, mat::parameters<col_major> > D(A);

    test(C, "Row-major matrix const");
    test(D, "Column-major matrix const");

    return 0;
}
 












