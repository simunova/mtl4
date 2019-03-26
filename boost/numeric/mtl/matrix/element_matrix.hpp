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

#ifndef MTL_ELEMENT_MATRIX_INCLUDE
#define MTL_ELEMENT_MATRIX_INCLUDE

namespace mtl { namespace mat {


template <typename Matrix, typename Rows, typename Cols>
struct element_matrix_t
{
    explicit element_matrix_t(const Matrix& matrix, const Rows& rows, const Cols& cols)
	: matrix(matrix), rows(rows), cols(cols)
    {}
    
    const Matrix&  matrix;
    const Rows&    rows;
    const Cols&    cols;
};
   

template <typename Matrix, typename Rows, typename Cols>
element_matrix_t<Matrix, Rows, Cols>
inline element_matrix(const Matrix& matrix, const Rows& rows, const Cols& cols)
{
    return element_matrix_t<Matrix, Rows, Cols>(matrix, rows, cols);
}

template <typename Matrix, typename Rows>
element_matrix_t<Matrix, Rows, Rows>
inline element_matrix(const Matrix& matrix, const Rows& rows)
{
    return element_matrix_t<Matrix, Rows, Rows>(matrix, rows, rows);
}


}} // namespace mtl::matrix

#endif // MTL_ELEMENT_MATRIX_INCLUDE
