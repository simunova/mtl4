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

#ifndef META_MATH_LEAST_SIGNIFICANT_ONE_BIT_INCLUDE
#define META_MATH_LEAST_SIGNIFICANT_ONE_BIT_INCLUDE

namespace meta_math {

template <unsigned long X>
struct least_significant_one_bit
{
    static const unsigned long value= ((X ^ (X-1)) + 1) >> 1;
};


} // namespace meta_math

#endif // META_MATH_LEAST_SIGNIFICANT_ONE_BIT_INCLUDE
