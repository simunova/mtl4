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

#ifndef MTL_VECTOR_IOTA_INCLUDE
#define MTL_VECTOR_IOTA_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {
  
/// Assigns sequentially increasing values to %vector v
template <typename Vector>
void iota(Vector& v, const typename Collection<Vector>::value_type offset= 0)
{
    vampir_trace<3013> tracer;
    for (typename Collection<Vector>::size_type i= 0; i < size(v); i++)
	v[i]= i + offset;
}

}} // namespace mtl::vector

#endif // MTL_VECTOR_IOTA_INCLUDE
