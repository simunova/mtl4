// Copyright 2006. Peter Gottschling, Matthias Troyer, Rolf Bonderer
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

#ifndef MATH_IS_INVERTIBLE_INCLUDE
#define MATH_IS_INVERTIBLE_INCLUDE

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

namespace math {

template <typename Operation, typename Element>
struct is_invertible_t
{
    // bool operator()(const Operation&, const Element&) const;
};


// By default all elements are invertible w.r.t. addition
// If only part of the elements are invertible it shall be handled by specialization of the type
// Whether invertibility is relevant at all shall be concrolled by the user with concept maps
template <typename Element>
struct is_invertible_t< add<Element>, Element >
  : public std::binary_function<add<Element>, Element, Element>
{
    bool operator() (const add<Element>&, const Element&) const 
    {
	return true;
    }
};


// By default all non-zero elements are invertible w.r.t. multiplication
// If another part of the elements or all elements are invertible it shall be handled by specialization of the type
// Whether invertibility is relevant at all shall be concrolled by the user with concept maps
template <typename Element>
struct is_invertible_t< mult<Element>, Element >
  : public std::binary_function<mult<Element>, Element, Element>
{
    bool operator() (const mult<Element>&, const Element& v) const 
    {
	return v != zero(v);
    }
};


// Function is shorter than typetrait-like functor
template <typename Operation, typename Element>
inline bool is_invertible(const Operation& op, const Element& v)
{
    return is_invertible_t<Operation, Element>() (op, v);
}


namespace detail {
    
    // Helper type whose operator returns true if v is not 0
    // 0 must be convertible into Element and Element must be EqualityComparable
    template <typename Operation, typename Element>
    struct non_zero_is_invertible_t
    {
	bool operator() (const Operation&, const Element& v)
	{
	    return v != Element(0);
	}
    };

} // namespace detail

} // namespace math

#endif // MATH_IS_INVERTIBLE_INCLUDE
