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

#ifndef MTL_STATIC_ASSERT_INCLUDE
#define MTL_STATIC_ASSERT_INCLUDE

#include <boost/static_assert.hpp>

namespace mtl {

#ifdef MTL_WITH_STATICASSERT
#  define MTL_STATIC_ASSERT(Condition, Message) static_assert(Condition, Message)
#else
#  define MTL_STATIC_ASSERT(Condition, Message) BOOST_STATIC_ASSERT(Condition)
#endif

} // namespace mtl

#endif // MTL_STATIC_ASSERT_INCLUDE
