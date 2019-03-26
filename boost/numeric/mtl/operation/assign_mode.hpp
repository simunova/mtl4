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

#ifndef MTL_ASSIGN_MODE_INCLUDE
#define MTL_ASSIGN_MODE_INCLUDE

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl { namespace assign {

/// Functor for assigning the result of a calculation (e.g. a sum) to a typically scalar value
struct assign_sum
{
    static const bool init_to_zero= true;

    template <typename T>
    static void init(T& v)
    {
	using math::zero;
	v= zero(v);
    }

    // Sets x with y when empty or not initialized
    template <typename T, typename U>
    static void set_empty(T& x, const U& y)
    {
	x= y;
    }

    // The first update sets the value and avoids the zeroing
    template <typename T, typename U>
    static void first_update(T& x, const U& y)
    {
	x= y;
    }

    template <typename T, typename U>
    static void update(T& x, const U& y)
    {
	x+= y;
    }

    // To be used like sfunctor
    template <typename T, typename U>
    static T& apply(T& x, const U& y)
    {
	return x= y;
    }

};

/// Functor for incrementing a typically scalar value with the result of a calculation (e.g. a sum)
struct plus_sum
{
    static const bool init_to_zero= false;

    template <typename T>
    static void init(T&) {}

    // Sets x with y when empty or not initialized
    template <typename T, typename U>
    static void set_empty(T& x, const U& y)
    {
	x= y;
    }

    template <typename T, typename U>
    static void first_update(T& x, const U& y)
    {
	x+= y;
    }

    template <typename T, typename U>
    static void update(T& x, const U& y)
    {
	x+= y;
    }

    // To be used like sfunctor
    template <typename T, typename U>
    static T& apply(T& x, const U& y)
    {
	return x+= y;
    }
};


/// Functor for incrementing a typically scalar value with the result of a calculation (e.g. a sum)
struct minus_sum
{
    static const bool init_to_zero= false;

    template <typename T>
    static void init(T&) {}

    // Sets x with y when empty or not initialized
    template <typename T, typename U>
    static void set_empty(T& x, const U& y)
    {
	x= -y;
    }

    template <typename T, typename U>
    static void first_update(T& x, const U& y)
    {
	x-= y;
    }

    template <typename T, typename U>
    static void update(T& x, const U& y)
    {
	x-= y;
    }

    // To be used like sfunctor
    template <typename T, typename U>
    static T& apply(T& x, const U& y)
    {
	return x-= y;
    }
};



/// Functor for multiplying a typically scalar value with the result of a calculation (e.g. a sum)
struct times_sum
{
    static const bool init_to_zero= false;

    template <typename T>
    static void init(T&) {}

    // Sets x with y when empty or not initialized
    template <typename T, typename U>
    static void set_empty(T& x, const U&)
    {
	x= T(0); 
    }

    template <typename T, typename U>
    static void first_update(T& x, const U& y)
    {
	x*= y;
    }

    template <typename T, typename U>
    static void update(T& x, const U& y)
    {
	x*= y;
    }

    // To be used like sfunctor
    template <typename T, typename U>
    static T& apply(T& x, const U& y)
    {
	return x*= y;
    }
};

/// Functor for dividing a typically scalar value with the result of a calculation (e.g. a sum)
struct divide_sum
{
    static const bool init_to_zero= false;

    template <typename T>
    static void init(T&) {}

    // Sets x with y when empty or not initialized
    template <typename T, typename U>
    static void set_empty(T& x, const U& y)
    {
        MTL_DEBUG_THROW_IF(y == U(0), division_by_zero());
	x= T(0); 
    }

    template <typename T, typename U>
    static void first_update(T& x, const U& y)
    {
        MTL_DEBUG_THROW_IF(y == U(0), division_by_zero());
	x/= y;
    }

    template <typename T, typename U>
    static void update(T& x, const U& y)
    {
        MTL_DEBUG_THROW_IF(y == U(0), division_by_zero());
	x/= y;
    }

    // To be used like sfunctor
    template <typename T, typename U>
    static T& apply(T& x, const U& y)
    {
	return x/= y;
    }
};

}} // namespace mtl::assign

#endif // MTL_ASSIGN_MODE_INCLUDE
