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

#ifndef META_MATH_LOG_2_INCLUDE
#define META_MATH_LOG_2_INCLUDE

#include <boost/numeric/meta_math/is_power_of_2.hpp>

namespace meta_math {

// Computes the logarithm to the basis 2
// Without testing if power of 2 it rounds values down to next integer
template <unsigned long X>
struct log_2
{
    // BOOST_STATIC_ASSERT(is_power_of_2_meta<X>::value);
    static const unsigned long tmp= X >> 1, value= log_2<tmp>::value + 1;
};

template <> struct log_2<1>
{
    static const unsigned long value= 0;
};

template <> struct log_2<0>
{
  // #error "Logarithm of 0 is undefined"
  BOOST_STATIC_ASSERT(true); // Logarithm of 0 is undefined
};


} // namespace meta_math

#endif // META_MATH_LOG_2_INCLUDE
