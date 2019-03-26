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

#ifndef MTL_WRAPPED_OBJECT_INCLUDE
#define MTL_WRAPPED_OBJECT_INCLUDE

namespace mtl {

/// Wrapper for objects
/** Can be referred in base classes to preserve initialization order
    and avoid ambiguities. Using the object as base class can cause
    name clashes and as member it cannot be initialized after all
    base classes with pedantic warnings. **/
template <typename View>
struct wrapped_object
{
    template <typename T>
    wrapped_object(T& x) : wrapped_object_member(x) {}

    template <typename T, typename U>
    wrapped_object(const T& x, const U& y) : wrapped_object_member(x, y) {}

    View wrapped_object_member; // ugly name avoid ambiguities in derived classes
};

} // namespace mtl

#endif // MTL_WRAPPED_OBJECT_INCLUDE
