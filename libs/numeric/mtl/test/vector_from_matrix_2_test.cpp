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

// #define MTL_VERBOSE_TEST  // to turn on the output

#include <iostream>
#include <cmath>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  
using namespace mtl;
   
template <typename Matrix>
void test(Matrix& A, const char* name)
{
    A= 1, 2, 3,
       4, 5, 6, 
       7, 8, 9;
    mtl::io::tout << "\n" << name << "  --- A = \n" << A; 

    typename mtl::ColumnInMatrix<Matrix>::type c0(A[iall][0]);
    c0[2]= 10.0;    
    MTL_THROW_IF(A[2][0] != 10.0, mtl::runtime_error("Matrix modification in column did not work"));

    typename mtl::ColumnInMatrix<Matrix>::type c1(clone(A[iall][1]));
    c1[2]= 11.0;
    MTL_THROW_IF(A[2][1] != 8.0, mtl::runtime_error("Matrix modied by a cloned column"));
   
    dense_vector<float> c2(clone(A[iall][2]));
    c2[2]= 11.5;
    MTL_THROW_IF(A[2][2] != 9.0, mtl::runtime_error("Matrix modied by a cloned column"));
 
    typename mtl::RowInMatrix<Matrix>::type r0(A[0][iall]);
    r0[2]= 12.0;    
    MTL_THROW_IF(A[0][2] != 12.0, mtl::runtime_error("Matrix modification in row did not work"));

    typename mtl::RowInMatrix<Matrix>::type r1(clone(A[1][iall]));
    r1[2]= 13.0;
    MTL_THROW_IF(A[1][2] != 6.0, mtl::runtime_error("Matrix modied by a cloned row"));
    
    dense_vector<float, mtl::vec::parameters<row_major> > r2(clone(A[2][iall]));
    r2[2]= 13.5;
    MTL_THROW_IF(A[2][2] != 9.0, mtl::runtime_error("Matrix modied by a cloned row"));

    mtl::io::tout << "A = \n" << A << std::endl; 
}
 


int main(int, char**)
{
    using namespace mtl;

    dense2D<float>                                 A(3, 3);
    dense2D<float, mat::parameters<col_major> > B(3, 3);

    test(A, "Row-major matrix   ");     
    test(B, "Column-major matrix");  

    return 0;
}
 












