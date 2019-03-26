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

#ifndef MTL_ARPREC_INCLUDE
#define MTL_ARPREC_INCLUDE

#ifdef MTL_HAS_ARPREC

#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>

#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/utility/true_copy.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

namespace mtl {

    /// Specialization for ARPREC complex numbers
    template <>
    struct Magnitude<mp_complex>
    {
	/// The associated type is ARPREC real
	typedef mp_real type;
    };

    namespace traits {

	template <>
	struct true_copy<mp_real_temp>
	{
	    typedef mp_real type;
	};
    }
} // namespace mtl

namespace math {

    template <>
    struct identity_t< add<mp_complex>, mp_complex > 
      : public std::binary_function<add<mp_complex>, mp_complex, mp_complex>
    { 
	mp_complex operator() (const add<mp_complex>&, const mp_complex& /*ref*/) const
	{
	    return mp_complex("0", "0");
	}
    };

} // math

#endif // MTL_HAS_ARPREC

#endif // MTL_ARPREC_INCLUDE
