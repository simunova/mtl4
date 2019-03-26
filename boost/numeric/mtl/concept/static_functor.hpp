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

#ifndef MTL_STATIC_FUNCTOR_INCLUDE
#define MTL_STATIC_FUNCTOR_INCLUDE

#include <complex>

namespace mtl {

#ifdef __GXX_CONCEPTS__

// Concept to specify static unary functors
auto concept StaticUnaryFunctor<typename T>
{
    typename argument_type = T::argument_type;
    typename result_type = T::result_type;

    static result_type apply(argument_type);
    result_type T::operator()(argument_type);
};

auto concept StaticBinaryFunctor<typename T>
{
    typename first_argument_type  = T::first_argument_type;
    typename second_argument_type = T::second_argument_type;
    typename result_type          = T::result_type;

    static result_type apply(first_argument_type, second_argument_type);
    result_type T::operator()(first_argument_type, second_argument_type);
};

#else  // now without concepts

/// Concept/Type-trait for static unary functors
/** This name is overloaded: when MTL4 is compiled with a concept-compiler
    StaticUnaryFunctor is a concept otherwise a type-trait.
**/
template <typename T>
struct StaticUnaryFunctor
{
    /// Associated type for argument, by default set to class's internal typedef (specialize if not present)
    typedef typename T::argument_type    argument_type;
    /// Associated type for result, by default set to class's internal typedef (specialize if not present)
    typedef typename T::result_type      result_type;
};

/// Concept/Type-trait for static binary functors
/** This name is overloaded: when MTL4 is compiled with a concept-compiler
    StaticBinaryFunctor is a concept otherwise a type-trait.
**/
template <typename T>
struct StaticBinaryFunctor
{
    /// Associated type for 1st argument, by default set to class's internal typedef (specialize if not present)
    typedef typename T::first_argument_type    first_argument_type;
    /// Associated type for 2nd argument, by default set to class's internal typedef (specialize if not present)
    typedef typename T::second_argument_type   second_argument_type;
    /// Associated type for result, by default set to class's internal typedef (specialize if not present)
    typedef typename T::result_type            result_type;
};

#endif  // __GXX_CONCEPTS__

} // namespace mtl

#endif // MTL_STATIC_FUNCTOR_INCLUDE
