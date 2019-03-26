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

using mtl::io::tout;

template <typename Matrix>
void access_test(const Matrix& A)
{
    MTL_THROW_IF(A[0][0] != 1, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(A[0][1] != 0, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(A[0][2] != 2, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(A[1][1] != 3, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(A[2][2] != 4, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(A[2][3] != 5, mtl::runtime_error("wrong value"));
}

template <typename Matrix>
void iteration_test(const Matrix& A)
{
  using mtl::tag::nz;
  namespace traits = mtl::traits;

  typename traits::row<Matrix>::type         row(A);
  typename traits::col<Matrix>::type         col(A);
  typename traits::const_value<Matrix>::type value(A); 
  typedef typename traits::range_generator<nz, Matrix>::type Cursor;

  for (Cursor cursor = mtl::begin<nz>(A), cend = mtl::end<nz>(A); cursor != cend; ++cursor) {
     tout << "A[" << row(*cursor) << "," << col(*cursor) << "] = " << value(*cursor) << std::endl;
     MTL_THROW_IF(row(*cursor) == 0 && col(*cursor) == 2 && value(*cursor) != 2.0, 
		  mtl::runtime_error("wrong value"))
  }

}

template <typename Matrix>
void copy_test(Matrix& A, unsigned dim1, unsigned dim2)
{
  Matrix B;
  laplacian_setup(B, dim1, dim2);
  
  typename Matrix::size_type  *row_ptr = B.address_major();
  typename Matrix::size_type  *col_ind = B.address_minor();
  typename Matrix::value_type *entries = B.address_data();

  A.change_dim(num_rows(B), num_cols(B));
  A = Matrix(num_rows(B), num_cols(B), B.nnz(), row_ptr, col_ind, entries);
  
  A *= 2.;
  MTL_THROW_IF(A[2][2] != 8, mtl::runtime_error("wrong value after scaling"));
  MTL_THROW_IF(B[2][2] != 4, mtl::runtime_error("wrong value after scaling, probably aliasing"));

  tout << "Matrix A:\n" << A;
  tout << "Matrix B:\n" << B;
}

template <typename Matrix>
void test(Matrix& A, const char* type)
{
    std::cout << "Matrix A [" << type << "]:\n" << A << std::endl;  
    access_test(A);
    iteration_test(A);
    copy_test(A, 3, 3);
}

int main(int, char**)
{
    using namespace mtl; using std::size_t;

    typedef compressed2D<double> Matrix;

    size_t dim1= 3, dim2= 4, nnz = 5;

    /*
     * Initialize Matrix
     *
     *     [1 0 2 0]
     * A = [0 3 0 0]
     *     [0 0 4 5]
     */
    size_t row_ptr[] = {0, 2, 3, 5}, col_ind[] = {0, 2, 1, 2, 3};
    double entries[] = {1., 2., 3., 4., 5.};

    Matrix A(dim1, dim2, nnz, row_ptr, col_ind, entries);
    test(A, "row-major");

    typedef mtl::mat::parameters<mtl::tag::col_major>  Parameters2;
    typedef compressed2D<double, Parameters2>             Matrix2;

    size_t col_ptr2[] = {0, 1, 2, 4, 5}, row_ind2[] = {0, 1, 0, 2, 2};
    double entries2[] = {1., 3., 2., 4., 5.};

    Matrix2 A2(dim1, dim2, nnz, col_ptr2, row_ind2, entries2);
    test(A2, "column-major");

    return 0;
}
