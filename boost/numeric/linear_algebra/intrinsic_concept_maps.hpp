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

#ifndef MATH_INTRINSIC_CONCEPT_MAPS_INCLUDE
#define MATH_INTRINSIC_CONCEPT_MAPS_INCLUDE

#include <concepts>

namespace math {


    // The following concepts are used to classify intrinsic arithmetic types.
    // The standard concepts already define the syntactic requirements,
    // i.e. the interface.
    // However they sey nothing about the semantics.
    // Therefore, user-defined types can model the syntactic/interface
    // requirements while still having a different mathematical behavior.
    // For that reason, we introduce concepts that are only used for intrinsic types.
    // For them we can define concept_maps regarding semantic behavior as monoids.

    concept IntrinsicSignedIntegral<typename T> 
      : std::SignedIntegralLike<T>
    {}

    concept IntrinsicUnsignedIntegral<typename T> 
      : std::UnsignedIntegralLike<T>
    {}

    concept IntrinsicFloatingPoint<typename T>
      : std::FloatingPointLike<T> 
    {}


    // Intrinsic types are chategorized:

    // concept_map IntrinsicSignedIntegral<char> {}; ???
    concept_map IntrinsicSignedIntegral<signed char> {};
    concept_map IntrinsicUnsignedIntegral<unsigned char> {};
    concept_map IntrinsicSignedIntegral<short> {};
    concept_map IntrinsicUnsignedIntegral<unsigned short> {};
    concept_map IntrinsicSignedIntegral<int> {};
    concept_map IntrinsicUnsignedIntegral<unsigned int> {};
    concept_map IntrinsicSignedIntegral<long> {};
    concept_map IntrinsicUnsignedIntegral<unsigned long> {};
    concept_map IntrinsicSignedIntegral<long long> {};
    concept_map IntrinsicUnsignedIntegral<unsigned long long> {};

    concept_map IntrinsicFloatingPoint<float> { }
    concept_map IntrinsicFloatingPoint<double> { }
        
    concept Integral<typename T> : std::CopyAssignable<T>     
    {
        requires std::HasPlus<T> && std::HasMinus<T> && std::HasMultiply<T> && std::HasDivide<T>
              && std::HasUnaryPlus<T> && std::HasNegate<T>;


        T& operator++(T&);
        T operator++(T& t, int) { T tmp(t); ++t; return tmp; }
        T& operator--(T&);
        T operator--(T& t, int) { T tmp(t); --t; return tmp; }

        requires std::Convertible<std::HasUnaryPlus<T>::result_type, T>
              && std::Convertible<std::HasNegate<T>::result_type, T>
              && std::Convertible<std::HasPlus<T>::result_type, T>
              && std::Convertible<std::HasMinus<T>::result_type, T>
              && std::Convertible<std::HasMultiply<T>::result_type, T>
              && std::Convertible<std::HasDivide<T>::result_type, T>;
        
        T& operator*=(T&, T);
        T& operator/=(T&, T);
        T& operator+=(T&, T);
        T& operator-=(T&, T);
    	
        requires std::HasComplement<T> && std::HasModulus<T> && std::HasBitAnd<T>
              && std::HasBitOr<T> && std::HasBitXor<T> && std::HasLeftShift<T> 
              && std::HasRightShift<T>;

        requires std::Convertible<std::HasComplement<T>::result_type, T>
              && std::Convertible<std::HasModulus<T>::result_type, T>
              && std::Convertible<std::HasBitAnd<T>::result_type, T>
              && std::Convertible<std::HasBitOr<T>::result_type, T>
              && std::Convertible<std::HasBitXor<T>::result_type, T>
              && std::Convertible<std::HasLeftShift<T>::result_type, T>
              && std::Convertible<std::HasRightShift<T>::result_type, T>;

        requires std::LessThanComparable<T> && std::EqualityComparable<T>;

        T& operator%=(T&, T);
        T& operator&=(T&, T);
        T& operator|=(T&, T);
        T& operator^=(T&, T);
        T& operator<<=(T&, T);
        T& operator>>=(T&, T);
    }


    template <typename T>
    requires IntrinsicUnsignedIntegral<T>
    concept_map Integral<T> {typedef T result_type;}

    template <IntrinsicSignedIntegral T>
    concept_map Integral<T> {typedef T result_type;}


} // namespace math

#endif // MATH_INTRINSIC_CONCEPT_MAPS_INCLUDE
