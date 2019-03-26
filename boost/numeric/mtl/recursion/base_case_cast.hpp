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

#ifndef MTL_BASE_CASE_CAST_INCLUDE
#define MTL_BASE_CASE_CAST_INCLUDE

#include <boost/numeric/mtl/recursion/base_case_matrix.hpp>
#include <boost/numeric/mtl/recursion/simplify_base_case_matrix.hpp>

namespace mtl { namespace recursion {


template <typename BaseCaseTest, typename Matrix>
typename base_case_matrix<Matrix, BaseCaseTest>::type inline
base_case_cast(Matrix const& matrix)
{
    return simplify_base_case_matrix(matrix, BaseCaseTest());
}


template <typename BaseCaseTest, typename Matrix>
typename base_case_matrix<Matrix, BaseCaseTest>::type inline
base_case_cast(Matrix& matrix)
{
    return simplify_base_case_matrix(matrix, BaseCaseTest());
}


}} // namespace mtl::recursion

#endif // MTL_BASE_CASE_CAST_INCLUDE
