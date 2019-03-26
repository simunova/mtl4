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

#ifndef MATH_POWER_INCLUDE
#define MATH_POWER_INCLUDE

#include <concepts>
#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <stdexcept>


namespace math {

    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	std::cout << "[Magma] ";
	if (n < 1) throw std::range_error("power [magma]: n must be > 0");

	Element value= a;
	for (; n > 1; --n)
	    value= op(value, a);
	return value;
    }

#if 0
    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires SemiGroup<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element multiply_and_square_horner(const Element& a, Exponent n, Op op) 
    {
	if (n < 1) throw std::range_error("mult&square Horner: n must be > 0");

        // Set mask to highest bit
        Exponent mask= 1 << (8 * sizeof(mask) - 1);

        // If this is a negative number right shift can insert 1s instead of 0s -> infinite loop
        // Therefore we take the 2nd-highest bit
        if (mask < 0)
	    mask= 1 << (8 * sizeof(mask) - 2);

        // Find highest 1 bit
        while(!bool(mask & n)) mask>>= 1;

        Element value= a;
        for (mask>>= 1; mask; mask>>= 1) {
	    value= op(value, value);
	    if (n & mask) 
		value= op(value, a);
        }
        return value;
    }

    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires SemiGroup<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	return multiply_and_square_horner(a, n, op);
    }
#endif


#if 1
    // With Horner scheme we can avoid recursion  
    // This one is more intuitive (I believe)      
    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires SemiGroup<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	std::cout << "[SemiGroup] ";
	if (n < 1) throw std::range_error("power [SemiGroup]: n must be > 0");

	Exponent half(n / 2);
        // If half is 0 then n must be 1 and the result is a
        if (half == 0)
	    return a;

        // Compute power of downward rounded exponent and "square" the result
        Element value= power(a, half, op);
        value= op(value, value);

        // If n is odd another operation with a is needed
        if (n & 1) 
	    value= op(value, a);
        return value;
    }
#endif



    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires Monoid<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element multiply_and_square(const Element& a, Exponent n, Op op) 
    {
	// Same as the simpler form except that the first multiplication is made before 
	// the loop and one squaring is saved this way
	if (n < 0) throw std::range_error("mult&square: n must be >= 0");
	
	using math::identity;
	Element value= bool(n & 1) ? Element(a) : Element(identity(op, a)), square= a;
	
	for (n>>= 1; n > 0; n>>= 1) {
	    square= op(square, square); 
	    if (n & 1) 
		value= op(value, square);
	}
	return value;  
    } 


    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires Monoid<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	std::cout << "[Monoid] ";
	return multiply_and_square(a, n, op);
    }




    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires PIMonoid<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	std::cout << "[PIMonoid] ";
	if (n < 0 && !is_invertible(op, a)) 
	    throw std::range_error("power [PIMonoid]: a must be invertible with n < 0");

	return n < 0 ? multiply_and_square(Element(inverse(op, a)), Exponent(-n), op)
	             : multiply_and_square(a, n, op);
    }

#if 1
    template <typename Op, std::Semiregular Element, Integral Exponent>
        requires Group<Op, Element> 
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	std::cout << "[Group] ";
	// For groups we don't need any range test

	return n < 0 ? multiply_and_square(Element(inverse(op, a)), Exponent(-n), op)
	             : multiply_and_square(a, n, op);
    }
#endif


#if 0
    template <typename Op, typename Element, typename Exponent>
        requires Group<Op, Element> 
              && Integral<Exponent>
              && std::Semiregular<Element>
              && std::Callable2<Op, Element, Element>
              && std::Convertible<std::Callable2<Op, Element, Element>::result_type, Element>
              && std::Semiregular<math::Inversion<Op, Element>::result_type>
              && std::HasNegate<Exponent>
              && math::Monoid<Op, math::Inversion<Op, Element>::result_type>
              && Integral< std::HasNegate<Exponent>::result_type>
              && std::Callable2<Op, math::Inversion<Op, Element>::result_type, 
				math::Inversion<Op, Element>::result_type>
              && std::Convertible<std::Callable2<Op, math::Inversion<Op, Element>::result_type, 
						 math::Inversion<Op, Element>::result_type>::result_type, 
				  math::Inversion<Op, Element>::result_type>
    inline Element power(const Element& a, Exponent n, Op op)
    {
	return n < 0 ? multiply_and_square(inverse(op, a), -n, op)
	             : multiply_and_square(a, n, op);
    }
#endif

} // namespace math

#endif // MATH_POWER_INCLUDE
