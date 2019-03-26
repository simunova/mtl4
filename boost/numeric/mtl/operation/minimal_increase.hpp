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

#ifndef MTL_MINIMAL_INCREASE_INCLUDE
#define MTL_MINIMAL_INCREASE_INCLUDE

#include <limits>

namespace mtl {

/// Increase x minimally: if x == 0 take minimal value, if x > 0 multiply by (1+2eps) otherwise divide by
template <typename T>
T inline minimal_increase(const T& x)
{
    const T factor= T(1) + T(2) * std::numeric_limits<T>::epsilon();
    if (x == T(0))
	return std::numeric_limits<T>::denorm_min();
    else 
	return x > T(0) ? x * factor : x / factor;	    
}



} // namespace mtl

#endif // MTL_MINIMAL_INCREASE_INCLUDE
