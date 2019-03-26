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

#ifndef MTL_MAGNITUDE_INCLUDE
#define MTL_MAGNITUDE_INCLUDE

#include <complex>

namespace mtl {

// =================================================
// Concept to specify return type of abs (and norms)
// =================================================


#ifdef __GXX_CONCEPTS__

// Concept to specify to specify projection of scalar value to comparable type
// For instance as return type of abs
// Minimalist definition for maximal applicability
auto concept Magnitude<typename T>
{
    typename type = T;
};

template <typename T>
concept_map Magnitude<std::complex<T> >
{
    typedef T type;
}


// Concept for norms etc., which are real values in mathematical definitions
auto concept RealMagnitude<typename T>
  : Magnitude<T>
{
    requires std::EqualityComparable<type>;
    requires std::LessThanComparable<type>;

#ifndef CONCEPTS_WITHOUT_OVERLOADED_REQUIREMENTS
    requires Field<type>;
#endif

    type sqrt(type);

    type abs(T);
}

#else  // now without concepts

/// Concept/Type-trait for magnitudes of scalar values
/** This name is overloaded: when MTL4 is compiled with a concept-compiler
    Magnitude is a concept otherwise a type-trait.
    It is used for instance in norms. 
**/
template <typename T>
struct Magnitude
{
    /// Associated type; the default is T; must be specialized appropriately
    typedef T type;
};

/// Specialization for complex numbers
template <typename T>
struct Magnitude<std::complex<T> >
{
    /// The associated type is defined to the complex's value type
    typedef T type;
};

template <typename T> struct RealMagnitude
  : public Magnitude<T>
{};

#endif  // __GXX_CONCEPTS__

} // namespace mtl

#endif // MTL_MAGNITUDE_INCLUDE
