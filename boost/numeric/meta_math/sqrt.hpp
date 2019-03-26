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

#ifndef META_MATH_SQRT_INCLUDE
#define META_MATH_SQRT_INCLUDE

#include <boost/numeric/meta_math/abs.hpp>

namespace meta_math {

template <long int root, long int x>
struct sqrt_check
{
    static bool const value = root * root <= x && (root+1) * (root+1) > x;
};


namespace impl {

    template <long int guess, long int x, bool Converged>
    struct sqrt_impl
    {
	typedef long int type;
	typedef sqrt_impl   self;
	static long int const quotient = x / guess,
	                      new_value = (quotient + guess) / 2;
	static bool const converging = abs<guess - quotient>::value < 2;
	static long int const value = sqrt_impl<new_value, x, converging>::value;
    };

    // If the condition becomes true the guessed root will be the returned value
    template <long int guess, long int x>
    struct sqrt_impl<guess, x, true> 
    {
	static long int const value = guess;
    };

}

template <long int x>
struct sqrt 
{
  typedef long int type;
  static long int const value = impl::sqrt_impl<1, x, false>::value;
};  


} // namespace meta_math

#endif // META_MATH_SQRT_INCLUDE
