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

#ifndef MTL_LAZY_INCLUDE
#define MTL_LAZY_INCLUDE

#include <boost/numeric/mtl/operation/assign_mode.hpp>

namespace mtl {

/// Helper class for lazy evaluation
template <typename T>
struct lazy_t
{
    lazy_t(T& data) : data(data) {}

    template <typename U>
    lazy_assign<T, U, assign::assign_sum> operator=(const U& other) 
    { return lazy_assign<T, U, assign::assign_sum>(data, other); }

    template <typename U>
    lazy_assign<T, U, assign::plus_sum> operator+=(const U& other) 
    { return lazy_assign<T, U, assign::plus_sum>(data, other); }

    template <typename U>
    lazy_assign<T, U, assign::minus_sum> operator-=(const U& other) 
    { return lazy_assign<T, U, assign::minus_sum>(data, other); }

    T& data;
};

/// Do not compute or assign x immediately but delay it until explicitly performed
template <typename T>
inline lazy_t<T> lazy(T& x) 
{ return lazy_t<T>(x); }

/// Do not compute or assign x immediately but delay it until explicitly performed
template <typename T>
inline lazy_t<const T> lazy(const T& x) 
{ return lazy_t<const T>(x); }



} // namespace mtl

#endif // MTL_LAZY_INCLUDE
