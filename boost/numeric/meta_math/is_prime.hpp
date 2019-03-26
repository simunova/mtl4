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

#ifndef META_MATH_IS_PRIME_INCLUDE
#define META_MATH_IS_PRIME_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/meta_math/sqrt.hpp>
// #include <boost/config/concept_macros.hpp>
#ifdef __GXX_CONCEPTS__
#   include <concepts>
#endif


namespace meta_math {

namespace mpl = boost::mpl;

namespace impl {

    // Checks if 'x' is divisible by some odd number >= max_odd, checking decreasingly
    template <long int x, long int max_odd>
    struct is_prime_to_max_odd
    {
	static bool const value = x % max_odd != 0
                                  && is_prime_to_max_odd<x, max_odd-2>::value;
    };

    // Once we reach 1, it's prime
    template <long int x> struct is_prime_to_max_odd<x, 1> : mpl::true_ {};


    // Returns the largest number that x is tried to divided by.
    // This is a odd number slightly larger than the approximated square root.
    template <long int x>
    struct max_prime_compare
    {
	static long int const tmp = sqrt<x>::value,
                              value = tmp % 2 == 0 ? tmp + 1 : tmp + 2;
    };


    // Checks if there is an odd number between 3 and sqrt(x)+1 that divides x
    // if there is no divisor found then x is prime (otherwise not)
    // must be disabled when x is even
    template <long int x, bool Disable>
    struct check_odd
    {
	static bool const value = is_prime_to_max_odd<x, max_prime_compare<x>::value>::value;
    };

    // Intended for even numbers (> 2) which are always prime
    template <long int x>
    struct check_odd<x, true>
    {
	static bool const value = false;
    };

}

template <long int x>
struct is_prime
{
  static bool const value = impl::check_odd<x, x % 2 == 0>::value;
};

template <> struct is_prime<0> : mpl::false_ {};
template <> struct is_prime<1> : mpl::false_ {};
template <> struct is_prime<2> : mpl::true_ {};
template <> struct is_prime<3> : mpl::true_ {};
template <> struct is_prime<5> : mpl::true_ {};


#ifdef __GXX_CONCEPTS__
    concept Prime<long int N> {}

    template <long int N>
        where std::True<is_prime<N>::value> 
    concept_map Prime<N> {}
#endif



} // namespace meta_math

#endif // META_MATH_IS_PRIME_INCLUDE
