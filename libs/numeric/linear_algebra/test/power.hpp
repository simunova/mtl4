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

#ifndef MTL_POWER_INCLUDE
#define MTL_POWER_INCLUDE

#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>
#include <iostream>

namespace mtl {

template <typename Op, typename Element, typename Exponent>
  _GLIBCXX_WHERE( math::Magma<Op, Element> 
            && std::Integral<Exponent> )      
inline Element power(const Element& a, Exponent n, Op op) 
{
#   ifdef MTL_TRACE_POWER_DISPATCHING 
       std::cout << "[Magma] ";
#   endif

    if (n < 1) throw "In power: exponent must be greater than 0";
    
    Element value= a;
    for (; n > 1; --n)
	value= op(value, a);
    return value;
}


# ifndef __GXX_CONCEPTS__
#   ifdef LA_SHOW_WARNINGS
#     warning "Automatic dispatching only works with concept compiler"
#     warning "If structure is a Monoid you can call square_and_multiply directly"
#   endif
# else

template <typename Op, typename Element, typename Exponent>
    where math::SemiGroup<Op, Element> && std::Integral<Exponent>
inline Element power(const Element& a, Exponent n, Op op)
{
#   ifdef MTL_TRACE_POWER_DISPATCHING 
       std::cout << "[SemiGroup] ";
#   endif

    return recursive_multiply_and_square(a, n, op);
}

template <typename Op, typename Element, typename Exponent>
    where math::Monoid<Op, Element> && std::Integral<Exponent>
inline Element power(const Element& a, Exponent n, Op op)
{
#   ifdef MTL_TRACE_POWER_DISPATCHING 
       std::cout << "[Monoid] ";
#   endif

    return multiply_and_square(a, n, op);
}

template <typename Op, typename Element, typename Exponent>
    where math::PartiallyInvertibleMonoid<Op, Element> && std::SignedIntegral<Exponent>
inline Element power(const Element& a, Exponent n, Op op)
{
#   ifdef MTL_TRACE_POWER_DISPATCHING 
       std::cout << "[PIMonoid] ";
#   endif
    using math::inverse; using math::is_invertible;

    if (n < 0 && !is_invertible(op, a)) 
        throw "In power [PIMonoid]: a must be invertible with negative n";

    return n >= 0 ? multiply_and_square(a, n, op) 
	          : multiply_and_square(inverse(op, a), -n, op);
}

template <typename Op, typename Element, typename Exponent>
    where math::Group<Op, Element> && std::SignedIntegral<Exponent>
inline Element power(const Element& a, Exponent n, Op op)
{
#   ifdef MTL_TRACE_POWER_DISPATCHING 
       std::cout << "[Group] ";
#   endif
    using math::inverse;

    return n >= 0 ? multiply_and_square(a, n, op) 
                  : multiply_and_square(inverse(op, a), -n, op);
}


# endif  // __GXX_CONCEPTS__

} // namespace mtl

#endif // MTL_POWER_INCLUDE
