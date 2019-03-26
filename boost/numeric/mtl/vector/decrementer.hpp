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

#ifndef MTL_VECTOR_DECREMENTER_INCLUDE
#define MTL_VECTOR_DECREMENTER_INCLUDE

namespace mtl { namespace vec {

struct decrementer_base {};

/// CRTP class to decrement a vector with the result of the derived class
template <typename Derived>
struct decrementer : decrementer_base
{
    /// Function that must be defined in \p Derived where it is called by static down-cast
    template <typename Vector>
    void decrement_it(Vector& tgt) const
    {
	static_cast<const Derived&>(*this).decrement_it(tgt);
    }

};

}} // namespace mtl::vector

#endif // MTL_VECTOR_DECREMENTER_INCLUDE
