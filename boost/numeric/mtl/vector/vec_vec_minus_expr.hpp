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


#ifndef MTL_VEC_VEC_MINUS_EXPR_INCLUDE
#define MTL_VEC_VEC_MINUS_EXPR_INCLUDE

#include <boost/static_assert.hpp>

#include <boost/numeric/mtl/vector/vec_vec_pmop_expr.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>

namespace mtl { namespace vec {

template <typename E1, typename E2>
inline vec_vec_pmop_expr< E1, E2, mtl::sfunctor::minus<typename E1::value_type, typename E2::value_type> >
operator- (const vec_expr<E1>& e1, const vec_expr<E2>& e2)
{
    // do not subtract row and column vectors (or inconsistent value types)
    MTL_STATIC_ASSERT((boost::is_same<typename ashape::ashape<E1>::type, 
			                typename ashape::ashape<E1>::type>::value),
		      "Vectors do not have consistent algebraic shape (i.e. nested types).");
    typedef vec_vec_pmop_expr< E1, E2, mtl::sfunctor::minus<typename E1::value_type, typename E2::value_type> > type;
    return type(static_cast<const E1&>(e1), static_cast<const E2&>(e2));
}

} } // Namespace mtl::vector

#endif

