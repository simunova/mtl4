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

#ifndef MTL_MTL_ASSERT_INCLUDE
#define MTL_MTL_ASSERT_INCLUDE

#include <cassert>

namespace mtl {

#define MTL_ASSERT(Cond, Msg) assert((Cond) && Msg);

#define MTL_CRASH_IF(Cond, Msg) assert(!(Cond) && Msg);

} // namespace mtl

#endif // MTL_MTL_ASSERT_INCLUDE
