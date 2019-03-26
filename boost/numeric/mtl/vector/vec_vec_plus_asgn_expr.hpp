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

// Adapted from GLAS implementation by Karl Meerbergen and Toon Knappen


#ifndef MTL_VEC_VEC_PLUS_ASGN_EXPR_INCLUDE
#define MTL_VEC_VEC_PLUS_ASGN_EXPR_INCLUDE

#include <boost/static_assert.hpp>

#include <boost/numeric/mtl/vector/vec_vec_aop_expr.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>

namespace mtl { namespace vec {

// Model of VectorExpression
template <class E1, class E2>
struct vec_vec_plus_asgn_expr 
  : public vec_vec_aop_expr< E1, E2, mtl::sfunctor::plus_assign<typename E1::value_type, typename E2::value_type> >
{
    typedef vec_vec_aop_expr< E1, E2, mtl::sfunctor::plus_assign<typename E1::value_type, typename E2::value_type> > base;
    vec_vec_plus_asgn_expr( E1& v1, E2 const& v2 )
      : base( v1, v2 )
    {}
};

} } // Namespace mtl::vector

#endif

