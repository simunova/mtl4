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

#ifndef MTL_TRAITS_ALGEBRAIC_CATEGORY_INCLUDE
#define MTL_TRAITS_ALGEBRAIC_CATEGORY_INCLUDE

#include <boost/numeric/mtl/utility/is_what.hpp>

namespace mtl { namespace traits {

/// Meta-function for categorizing types into tag::scalar, tag::vector, and tag::matrix
/** Automatically derived from category 
    @ingroup Tags
*/
template <typename T>
struct algebraic_category
  : boost::mpl::if_<
	is_matrix<T>
      , tag::matrix
      , typename boost::mpl::if_<
       	    is_vector<T>
	  , tag::vector
	  , tag::scalar
	>::type
    >
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_ALGEBRAIC_CATEGORY_INCLUDE
