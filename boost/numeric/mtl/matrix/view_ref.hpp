// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_MATRIX_VIEW_REF_INCLUDE
#define MTL_MATRIX_VIEW_REF_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/hermitian_view.hpp>


namespace mtl { namespace mat {

template <typename Matrix>
inline Matrix& view_ref(Matrix& A)
{    return A; }

template <typename Matrix>
inline const Matrix& view_ref(const Matrix& A)
{    return A; }

template <typename Matrix>
inline Matrix& view_ref(transposed_view<Matrix>& A)
{    return A.ref; }

template <typename Matrix>
inline const Matrix& view_ref(transposed_view<const Matrix>& A)
{    return A.ref; }

template <typename Matrix>
inline const Matrix& view_ref(const transposed_view<Matrix>& A)
{    return A.ref; }

template <typename Matrix>
inline const Matrix& view_ref(const conj_view<Matrix>& A)
{    return A.ref; }

template <typename Matrix>
inline const Matrix& view_ref(const hermitian_view<Matrix>& A)
{    return A.const_ref(); }


}} // namespace mtl::matrix

#endif // MTL_MATRIX_VIEW_REF_INCLUDE
