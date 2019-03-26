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

#ifndef MTL_COPYSIGN_INCLUDE
#define MTL_COPYSIGN_INCLUDE

#include <complex>
#include <cmath>
#include <boost/numeric/mtl/operation/real.hpp>
#include <boost/numeric/mtl/operation/imag.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

namespace sfunctor {

    template <typename Value1, typename Value2>
    struct copysign
    {
	static inline Value1 apply(const Value1& v, const Value2& s)
	{
	    using math::zero; using std::abs;

	    Value1 a(abs(v));
	    return s < zero(s) ? -a : a;
	}
    };

    // This specialization is questionable and thus subject to elimination 
    template <typename Value1, typename Value2>
    struct copysign<std::complex<Value1>, Value2>
    {
	static inline Value1 apply(const std::complex<Value1>& v, const Value2& s)
	{
	    using mtl::real; using mtl::imag;
	    return std::complex<Value1>(copysign<Value1, Value2>::apply(real(v), s),
					copysign<Value1, Value2>::apply(imag(v), s));
	}
    };
#if 0 // ndef _MSC_VER doesn't work on BigRed with gcc 4.2.2 (????) and isn't used anyway
    template <>
    struct copysign<float, float>
    {
	static inline float apply(float v, float s)
	{ return ::copysignf(v, s); }
    };

    template <>
    struct copysign<double, double>
    {
	static inline double apply(double v, double s)
	{ return ::copysign(v, s); }
    };

    template <>
    struct copysign<long double, long double>
    {
	static inline long double apply(long double v, long double s)
	{ return copysignl(v, s); }
    };
#endif // _MSC_VER
}

/// sign of scalars; for complex numbers sign of real part
template <typename Value1, typename Value2>
inline Value1 copysign(const Value1& v, const Value2& s)
{
    vampir_trace<1001> tracer;
    return sfunctor::copysign<Value1, Value2>::apply(v, s);
}




} // namespace mtl

#endif // MTL_COPYSIGN_INCLUDE
