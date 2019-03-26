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

#ifndef MTL_POWER_EXPLICIT_INCLUDE
#define MTL_POWER_EXPLICIT_INCLUDE

#include <boost/numeric/linear_algebra/algebraic_concepts.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>



namespace mtl {

template <typename Op, typename Element, typename Exponent>
  _GLIBCXX_WHERE( std::Integral<Exponent> 
	    && std::Callable2<Op, Element, Element>
	    && std::Assignable<Element, std::Callable2<Op, Element, Element>::result_type>)            
inline Element power(const Element& base, Exponent n, Op op) 
{
    if (n < 1) throw "In power: exponent must be greater than 0";
    // std::cout << "[Magma] ";
    
    Element value= base;
    for (; n > 1; --n)
	value= op(value, base);
    return value;
}


# ifndef __GXX_CONCEPTS__
#   ifdef LA_SHOW_WARNINGS
#     warning "Automatic dispatching only works with concept compiler"
#     warning "If structure is a Monoid you can call square_and_multiply directly"
#   endif
# else

template <typename Op, typename Element, typename Exponent>
    where algebra::SemiGroup<Op, Element> && std::Integral<Exponent>
          && std::Callable2<Op, Element, Element>
          && std::Assignable<Element, std::Callable2<Op, Element, Element>::result_type>            
inline Element power(const Element& base, Exponent n, Op op)
{
    // std::cout << "[SemiGroup] ";

    if (n <= 0) throw "In recursive_multiply_and_square: exponent must greater than 0";

    Exponent half= n >> 1;

    // If halt is 0 then n must be 1 and the result is base
    if (half == 0)
	return base;

    // compute power of downward rounded exponent and square the result
    Element value= power(base, half, op);
    value= op(value, value);

    // if odd another multiplication with base is needed
    if (n & 1) 
	value= op(value, base);
    return value;
}

// {Op, Element} must be a Monoid
template <typename Op, typename Element, typename Exponent>
    where algebra::Monoid<Op, Element> 
          && std::Integral<Exponent>
          && std::Callable2<Op, Element, Element>
          && std::Assignable<Element, std::Callable2<Op, Element, Element>::result_type>
// && std::Assignable<Element, algebdra::Monoid<Op, Element>::identity_result_type>
          && std::Assignable<Element, Element>
inline Element multiply_and_square(const Element& base, Exponent n, Op op) 
{
    // Same as the simpler form except that the first multiplication is made before 
    // the loop and one squaring is saved this way
    if (n < 0) throw "In multiply_and_square: negative exponent";

    using math::identity;
    Element value= identity(op, base), square= identity(op, base);

    if (n & 1)
        value= base;

    for (n>>= 1; n > 0; n>>= 1) {
	square= op(square, square); 
	if (n & 1) 
	    value= op(value, square);
    }
    return value;  
} 

template <typename Op, typename Element, typename Exponent>
    where algebra::Monoid<Op, Element> && std::Integral<Exponent>
          && std::Callable2<Op, Element, Element>
          && std::Assignable<Element, std::Callable2<Op, Element, Element>::result_type>            
          && std::Assignable<Element, Element>
inline Element power(const Element& base, Exponent n, Op op)
{
    return multiply_and_square(base, n, op);
}

template <typename Op, typename Element, typename Exponent>
    where algebra::Group<Op, Element> && std::SignedIntegral<Exponent>
          && std::Callable2<Op, Element, Element>
          && std::Assignable<Element, std::Callable2<Op, Element, Element>::result_type>            
          && std::Assignable<Element, Element>
inline Element power(const Element& base, Exponent n, Op op)
{
    using math::inverse;

    return n >= 0 ? multiply_and_square(base, n, op) 
	          : multiply_and_square(inverse(op, base), -n, op);
}


# endif   // __GXX_CONCEPTS__

} // namespace mtl

#endif // MTL_POWER_EXPLICIT_INCLUDE
