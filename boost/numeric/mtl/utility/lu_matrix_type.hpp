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

#ifndef MTL_TRAITS_LU_MATRIX_TYPE_INCLUDE
#define MTL_TRAITS_LU_MATRIX_TYPE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

template <typename Matrix>
struct lu_matrix_type
{
    typedef Matrix type;
};

template <typename Value, typename Parameters> 
struct lu_matrix_type<mat::compressed2D<Value, Parameters> >
{
    typedef mat::dense2D<Value, Parameters> type;
};


}} // namespace mtl::traits

#endif // MTL_TRAITS_LU_MATRIX_TYPE_INCLUDE
