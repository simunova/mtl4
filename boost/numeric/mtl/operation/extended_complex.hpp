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

#ifndef MTL_OPERATION_EXTENDED_COMPLEX_INCLUDE
#define MTL_OPERATION_EXTENDED_COMPLEX_INCLUDE

#include <complex>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/numeric/mtl/utility/different_non_complex.hpp>
#include <boost/numeric/mtl/utility/extended_complex.hpp>


// Now we dare writing in the standard namespace
namespace std {

// ========
// Addition
// ========

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator+(const T& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(x + real(y), imag(y));
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator+(const complex<T>& x, const U& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) + y, imag(x));
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator+(const complex<T>& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) + real(y), imag(x) + imag(y));
}

// ===========
// Subtraction
// ===========

namespace detail {
    
    // To avoid stupid warnings on unary - for unsigned int specialize it
    template <typename T>
    typename boost::disable_if<boost::is_unsigned<T>, T>::type
    inline negate_helper(const T& x) 
    { return -x; }

    template <typename T>
    typename boost::enable_if<boost::is_unsigned<T>, T>::type
    inline negate_helper(const T& x) 
    { return T(0) - x; }
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator-(const T& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(x - real(y), detail::negate_helper(imag(y))); // see above for helper
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator-(const complex<T>& x, const U& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) - y, imag(x));
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator-(const complex<T>& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) - real(y), imag(x) - imag(y));
}

// ==============
// Multiplication
// ==============

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator*(const T& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(x * real(y), x * imag(y));
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator*(const complex<T>& x, const U& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) * y, imag(x) * y);
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator*(const complex<T>& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) * real(y) - imag(x) * imag(y),
		real(x) * imag(y) + imag(x) * real(y));
}

// ========
// Division
// ========
// Not necessarily the most efficient implementations
// Dealing primarily with types

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator/(const T& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    type r(x);
    return r/= type(real(y), imag(y));    
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator/(const complex<T>& x, const U& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    return type(real(x) / y, imag(x) / y);
}

template <typename T, typename U>
typename mtl::traits::extended_complex<T, U>::type // implicit enable_if
inline operator/(const complex<T>& x, const complex<U>& y)
{
    typedef typename mtl::traits::extended_complex<T, U>::type type;
    type r(real(x), imag(x));
    return r/= type(real(y), imag(y));    
}

     
// ==========
// Comparison
// ==========

template <typename T, typename U>
typename boost::enable_if<mtl::traits::different_non_complex<T, U>, bool>::type
inline operator==(const complex<T>& x, const U& y)
{
    return real(x) == y && imag(x) == T(0);
}

template <typename T, typename U>
typename boost::enable_if<mtl::traits::different_non_complex<T, U>, bool>::type
inline operator==(const T& x, const complex<U>& y)
{
    return x == real(y) && imag(y) == T(0);
}

template <typename T, typename U>
typename boost::enable_if<mtl::traits::different_non_complex<T, U>, bool>::type
inline operator==(const complex<T>& x, const complex<U>& y)
{
    return real(x) == real(y) && imag(x) == imag(y);
}

template <typename T, typename U>
typename boost::enable_if<mtl::traits::different_non_complex<T, U>, bool>::type
inline operator!=(const complex<T>& x, const U& y)
{
    return real(x) != y || imag(x) != T(0);
}

template <typename T, typename U>
typename boost::enable_if<mtl::traits::different_non_complex<T, U>, bool>::type
inline operator!=(const T& x, const complex<U>& y)
{
    return x != real(y) || imag(y) != T(0);
}

template <typename T, typename U>
typename boost::enable_if<mtl::traits::different_non_complex<T, U>, bool>::type
inline operator!=(const complex<T>& x, const complex<U>& y)
{
    return real(x) != real(y) || imag(x) != imag(y);
}

} // namespace std

#endif // MTL_OPERATION_EXTENDED_COMPLEX_INCLUDE
