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

#ifndef META_MATH_IS_POWER_OF_2_INCLUDE
#define META_MATH_IS_POWER_OF_2_INCLUDE

#include <boost/numeric/meta_math/least_significant_one_bit.hpp>

namespace meta_math {

template <unsigned long X>
struct is_power_of_2
{
    static const bool value= X == least_significant_one_bit<X>::value;
};

} // namespace meta_math

#endif // META_MATH_IS_POWER_OF_2_INCLUDE
