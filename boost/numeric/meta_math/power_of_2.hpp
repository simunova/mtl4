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

#ifndef META_MATH_POWER_OF_2_INCLUDE
#define META_MATH_POWER_OF_2_INCLUDE

namespace meta_math {


// Computes the n-th power of 2
// So simple, everybody could do it, it is only there for the sake of completeness
template <unsigned long N>
struct power_of_2
{
    static const unsigned long value= 1 << N;
};


} // namespace meta_math

#endif // META_MATH_POWER_OF_2_INCLUDE
