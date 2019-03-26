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

#ifndef MTL_MATRIX_SPLIT_COMPLEX_MATRIX_INCLUDE
#define MTL_MATRIX_SPLIT_COMPLEX_MATRIX_INCLUDE

#if 0 // Do we really need this???

namespace mtl { namespace matrix {

// Split one complex-valued matrix into two real-valued matrices.
/* Elements of the real matrix must be assignable from the real and imaginary part of the complex elements.
    Real matrices are resized if their size is 0 otherwise the matrices must have
    the same dimension. **/
template <typename MatrixComplex, typename MatrixReal, typename MatrixImaginary>
inline void split_complex_matrix(const MatrixComplex& c, MatrixReal& r, MatrixImaginary& i)
{
    r.checked_change_dim(num_rows(c), num_cols(c));
    i.checked_change_dim(num_rows(c), num_cols(c));

    c= 
    

}


}} // namespace mtl::matrix

#endif

#endif // MTL_MATRIX_SPLIT_COMPLEX_MATRIX_INCLUDE
