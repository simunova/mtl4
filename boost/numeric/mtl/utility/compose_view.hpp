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

#ifndef MTL_TRAITS_COMPOSE_VIEW_INCLUDE
#define MTL_TRAITS_COMPOSE_VIEW_INCLUDE

#include <boost/type_traits/remove_reference.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>

namespace mtl { namespace traits {

/// Compose matrix view from code and matrix type
template <unsigned Code, typename Matrix>
struct matrix_compose_view
{
    MTL_STATIC_ASSERT(Code <= 7, "Illegal ViewCode!");
    typedef typename boost::remove_reference<typename matrix_compose_view<Code & 1, Matrix>::type>::type matrix_type;

    typedef typename boost::mpl::if_c<
	(Code >= 6),
	mtl::mat::hermitian_view<matrix_type>,
	typename boost::mpl::if_c<
	    (Code >= 4),
	    mtl::mat::transposed_view<matrix_type>,
	    mtl::mat::conj_view<matrix_type>
	   >::type
       >::type type;
};

template <typename Matrix>
struct matrix_compose_view<0, Matrix>
{
    typedef Matrix& type;
};

template <typename Matrix>
struct matrix_compose_view<1, Matrix>
{
    typedef const Matrix& type;
};

// Add vector stuff

/// Compose view from code and collection type (matrix or vector)
template <unsigned Code, typename Collection>
struct compose_view
  : matrix_compose_view<Code, Collection>
{
    MTL_STATIC_ASSERT(is_matrix<Collection>::value, "Currently only matrices are supported.");
};

}} // namespace mtl::traits

#endif // MTL_TRAITS_COMPOSE_VIEW_INCLUDE
