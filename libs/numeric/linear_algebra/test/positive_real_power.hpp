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

#ifndef MTL_POSITIVE_REAL_POWER_INCLUDE
#define MTL_POSITIVE_REAL_POWER_INCLUDE
 
#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>

#include <libs/numeric/linear_algebra/test/positive_real.hpp>

// User defined data types and operators

using mtl::positive_real;

struct magma_mult : public math::mult<positive_real> {};
struct semigroup_mult : public math::mult<positive_real> {};
struct monoid_mult : public math::mult<positive_real> {};
struct pim_mult : public math::mult<positive_real> {};     // Partially invertible monoid
struct group_mult : public math::mult<positive_real> {};

namespace math {
    template<> struct identity_t<monoid_mult, positive_real> 
	: public identity_t<mult<positive_real>, positive_real> {};

    template<> struct identity_t<pim_mult, positive_real> 
	: public identity_t<mult<positive_real>, positive_real> {};
    template<> struct inverse_t<pim_mult, positive_real> 
	: public inverse_t<mult<positive_real>, positive_real> {};
    template<> struct is_invertible_t<pim_mult, positive_real> 
	: public detail::non_zero_is_invertible_t<pim_mult, positive_real> {};

    template<> struct identity_t<group_mult, positive_real> 
	: public identity_t<mult<positive_real>, positive_real> {};
    template<> struct inverse_t<group_mult, positive_real> 
	: public inverse_t<mult<positive_real>, positive_real> {};
    template<> struct is_invertible_t<group_mult, positive_real> 
	: public is_invertible_t<mult<positive_real>, positive_real> {};
}


template <typename Op>
void compute_power(positive_real base, int exp, Op op, const char* structure)
{
    using mtl::power;
    try {
	std::cout << base << "^" << exp << " as " << structure << "  " << power(base, exp, op) << '\n';
    } catch (char const* message) {
	std::cout << "\n== Exception: " << message << '\n';
    }
}



#endif // MTL_POSITIVE_REAL_POWER_INCLUDE
