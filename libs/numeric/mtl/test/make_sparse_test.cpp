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
 
#define MTL_VERBOSE_TEST

#include <boost/numeric/mtl/mtl.hpp>

// Commands in Matlab
// S = sparse(i,j,s,m,n,nzmax)
// S = sparse(i,j,s,m,n)
// S = sparse(i,j,s)
// S = sparse(m,n)

using mtl::io::tout;

typedef mtl::compressed2D<double, mtl::mat::unsigned_parameters>         matrix_type;

template <typename Matrix>
inline void check(const Matrix& A, unsigned m, unsigned n)
{
    MTL_THROW_IF(num_rows(A) != m, mtl::unexpected_result("Not exact number of rows."));
    MTL_THROW_IF(num_cols(A) != n, mtl::unexpected_result("Not exact number of columns."));
    MTL_THROW_IF(A.nnz() != 5, mtl::unexpected_result("Not exact number of non-zeros."));

    MTL_THROW_IF(A[4][3] != 3.0, mtl::unexpected_result("A[4][3] should be 3."));
    MTL_THROW_IF(A[1][6] != 4.0, mtl::unexpected_result("A[1][6] should be 4."));
    MTL_THROW_IF(A[1][5] != 0.0, mtl::unexpected_result("A[1][5] should be empty."));
}

template <typename SV, typename VV>
inline void test1(const SV& rows, const SV& cols, const VV& values, unsigned m, unsigned n)
{
    matrix_type A;
    A= make_sparse(rows, cols, values, m, n);
    tout << "test1: A is\n" << A;

    check(A, m, n);
}

template <typename SV, typename VV>
inline void test2(const SV& rows, const SV& cols, const VV& values)
{
    matrix_type A;
    A= make_sparse(rows, cols, values);
    tout << "test2: A is\n" << A;

    check(A, 5, 7);
}

inline void test3(unsigned m, unsigned n)
{
    matrix_type A;
    A= mtl::make_sparse(m, n);
    tout << "test3: A is\n" << A;

    MTL_THROW_IF(num_rows(A) != m, mtl::unexpected_result("Not exact number of rows."));
    MTL_THROW_IF(num_cols(A) != n, mtl::unexpected_result("Not exact number of columns."));
    MTL_THROW_IF(A.nnz() != 0, mtl::unexpected_result("Not exact number of non-zeros."));

    MTL_THROW_IF(A[1][5] != 0.0, mtl::unexpected_result("A[1][5] should be empty."));
}

int main(int, char**)
{
    mtl::dense_vector<unsigned>     rows(7), cols(7);
    rows= 2, 0, 4, 1, 3, 1, 0;
    cols= 0, 3, 3, 6, 5, 6, 2;

    mtl::dense_vector<double>       values(7);
    values= 1, 2, 3, 3, 5, 1, 0;

    test1(rows, cols, values, 6, 8);
    test2(rows, cols, values);
    test3(6, 8);
    
    return 0;
}
