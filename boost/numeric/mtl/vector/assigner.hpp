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

#ifndef MTL_VECTOR_ASSIGNER_INCLUDE
#define MTL_VECTOR_ASSIGNER_INCLUDE

namespace mtl { namespace vec {

struct assigner_base {};

/// CRTP class to assign the result of the derived class to a vector
template <typename Derived>
struct assigner : assigner_base
{
    /// Function that must be defined in \p Derived where it is called by static down-cast
    template <typename Vector>
    void assign_to(Vector& tgt) const
    {
	static_cast<const Derived&>(*this).assign_to(tgt);
    }

};

}} // namespace mtl::vector

#endif // MTL_VECTOR_ASSIGNER_INCLUDE
