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

#ifndef MTL_SQUARED_ABS_INCLUDE
#define MTL_SQUARED_ABS_INCLUDE

#include <cmath>
#include <complex>
#include <boost/mpl/or.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

/// When squaring magnitudes of intrinsic types the abs can be omitted 
template <typename T>
typename boost::enable_if<boost::mpl::or_<boost::is_integral<T>, boost::is_floating_point<T> >, T>::type
inline squared_abs(const T& x)
{
    vampir_trace<10> tracer;
    return x * x;
}

/// Squaring complex numbers can be computed without square root
template <typename T>
T inline squared_abs(const std::complex<T>& z)
{
    vampir_trace<11> tracer;
    using std::norm;
    return norm(z);
}

/// When squaring magnitudes of non-intrinsic and non-complex types we use abs for security 
template <typename T>
typename boost::disable_if<boost::mpl::or_<boost::is_integral<T>, boost::is_floating_point<T> >, T>::type
inline squared_abs(const T& x)
{
    vampir_trace<12> tracer;
    using std::abs;
    T a= abs(x);
    return a * a;
}


} // namespace mtl

#endif // MTL_SQUARED_ABS_INCLUDE
