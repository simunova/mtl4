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

#ifndef LA_CONCEPTS_INCLUDE
#define LA_CONCEPTS_INCLUDE

#include <boost/config/concept_macros.hpp>

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else
#  ifdef LA_SHOW_WARNINGS
#    warning "Concepts are not used"
#  endif
#endif


#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/is_invertible.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/algebraic_concepts.hpp>
#include <complex>

// If desired one can disable the default concept maps with LA_NO_CONCEPT_MAPS

// We consider to change the namespace from math to numeric
// More precisely, the concepts may be moved into namespace numeric and the standard functions stay in math

/// Namespace for mathematical concepts
/** In contrast to the ones in algebra the concepts can require 
    basic implementation concepts like std::CopyAssignable */  
namespace math {

#ifdef __GXX_CONCEPTS__

// ==================================
// Classification of Arithmetic Types
// ==================================

// We omit this now because it became part of the standard concepts

#if 0

// In addtion to std::Integral
concept Float<typename T> 
  : std::DefaultConstructible<T>, std::CopyConstructible<T>,
    std::LessThanComparable<T>, std::EqualityComparable<T>
{
  T operator+(T);
  T operator+(T, T);
  T& operator+=(T&, T);
  T operator-(T, T);
  T operator-(T);
  T& operator-=(T&, T);
  T operator*(T, T);
  T& operator*=(T&, T);
  T operator/(T, T);
  T& operator/=(T&, T);

  requires std::CopyAssignable<T>
        && std::SameType<std::CopyAssignable<T>::result_type, T&>;
}

concept_map Float<float> {}
concept_map Float<double> {}
concept_map Float<long double> {}

// The difference to Float is the lack of LessThanComparable
concept Complex<typename T> 
  : std::DefaultConstructible<T>, std::CopyConstructible<T>,
    std::EqualityComparable<T>
{
  T operator+(T);
  T operator+(T, T);
  T& operator+=(T&, T);
  T operator-(T, T);
  T operator-(T);
  T& operator-=(T&, T);
  T operator*(T, T);
  T& operator*=(T&, T);
  T operator/(T, T);
  T& operator/=(T&, T);

  requires std::CopyAssignable<T> 
        && std::SameType<std::CopyAssignable<T>::result_type, T&>;
}

template <typename T>
  requires Float<T>
concept_map Complex<std::complex<T> > {}


// TBD: Concept Arithmetic is useless like this, it should have operations and then be the base for
// Integral, Float and Complex
concept Arithmetic<typename T> {}

template <typename T>
  requires std::Integral<T>
concept_map Arithmetic<T> {}

template <typename T>
  requires Float<T>
concept_map Arithmetic<T> {}

template <typename T>
  requires Arithmetic<T>
concept_map Arithmetic< std::complex<T> > {}

#endif 

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

concept_map IntrinsicSignedIntegral<char> {} 
concept_map IntrinsicSignedIntegral<signed char> {}
concept_map IntrinsicUnsignedIntegral<unsigned char> {}
concept_map IntrinsicSignedIntegral<short> {}
concept_map IntrinsicUnsignedIntegral<unsigned short> {}
concept_map IntrinsicSignedIntegral<int> {}
concept_map IntrinsicUnsignedIntegral<unsigned int> {}
concept_map IntrinsicSignedIntegral<long> {}
concept_map IntrinsicUnsignedIntegral<unsigned long> {}
concept_map IntrinsicSignedIntegral<long long> {}
concept_map IntrinsicUnsignedIntegral<unsigned long long> {}

concept_map IntrinsicFloatingPoint<float> {}
concept_map IntrinsicFloatingPoint<double> {}



// ================
// Utility Concepts
// ================

// Concepts for functions mapping to same type or convertible
auto concept UnaryIsoFunction<typename Operation, typename Element>
{
    requires std::Callable1<Operation, Element>;
    requires std::Convertible<std::Callable1<Operation, Element>::result_type, Element>;

    typename result_type = std::Callable1<Operation, Element>::result_type;
};


auto concept BinaryIsoFunction<typename Operation, typename Element>
{
    requires std::Callable2<Operation, Element, Element>;
    requires std::Convertible<std::Callable2<Operation, Element, Element>::result_type, Element>;

    typename result_type = std::Callable2<Operation, Element, Element>::result_type;
};

#if 0
auto concept CompatibleBinaryFunction<typename A1, typename A2, typename Result>
{
    typename result_type;
    result_type F(A1, A2);
    requires std::Convertible<result_type, Result>;
}
#endif

// ==================
// Algebraic Concepts
// ==================


auto concept Magma<typename Operation, typename Element>
    : BinaryIsoFunction<Operation, Element>
{
    requires std::CopyAssignable<Element>
	  && std::CopyAssignable<Element, BinaryIsoFunction<Operation, Element>::result_type>;
};


// For algebraic structures that are commutative but not associative
// As an example floating point numbers are commutative but not associative
//   w.r.t. addition and multiplication
auto concept CommutativeMagma<typename Operation, typename Element>
  : Magma<Operation, Element>, 
    algebra::Commutative<Operation, Element>
{};


// SemiGroup is a refinement which must be nominal
auto concept SemiGroup<typename Operation, typename Element>
  : Magma<Operation, Element>, 
    algebra::SemiGroup<Operation, Element>
{};


auto concept CommutativeSemiGroup<typename Operation, typename Element>
  : SemiGroup<Operation, Element>,
    CommutativeMagma<Operation, Element>
{};


// Adding identity
// auto 
concept Monoid<typename Operation, typename Element>
  : SemiGroup<Operation, Element>, 
    algebra::Monoid<Operation, Element> 
{
    requires std::Convertible<identity_result_type, Element>;
};


auto concept CommutativeMonoid<typename Operation, typename Element>
  : CommutativeSemiGroup<Operation, Element>, 
    Monoid<Operation, Element>
{};


concept PartiallyInvertibleMonoid<typename Operation, typename Element>
  : Monoid<Operation, Element>, 
    algebra::Inversion<Operation, Element> 
{
    typename is_invertible_result_type;
    is_invertible_result_type is_invertible(Operation, Element);
    requires std::Convertible<is_invertible_result_type, bool>;

    requires std::Convertible<inverse_result_type, Element>;

    // Does it overwrites the axiom from algebra::Inversion
    axiom Inversivity(Operation op, Element x)
    {
	// Only for invertible elements:
	if (is_invertible(op, x))
	    op( x, inverse(op, x) ) == identity(op, x); 
	if ( is_invertible(op, x) )
	    op( inverse(op, x), x ) == identity(op, x); 
    }
};


auto concept PartiallyInvertibleCommutativeMonoid<typename Operation, typename Element>
  : PartiallyInvertibleMonoid<Operation, Element>, 
    CommutativeMonoid<Operation, Element>   
{};


concept Group<typename Operation, typename Element>
  : PartiallyInvertibleMonoid<Operation, Element>,
    algebra::Group<Operation, Element>
{
    axiom AlwaysInvertible(Operation op, Element x)
    {
	is_invertible(op, x);
    }

    axiom GlobalInversivity(Operation op, Element x)
    {
	// In fact this is implied by AlwaysInvertible and inherited Inversion axiom
	// However, we don't rely on the compiler to deduce this
	op( x, inverse(op, x) ) == identity(op, x);
	op( inverse(op, x), x ) == identity(op, x);
    }
};


auto concept AbelianGroup<typename Operation, typename Element>
  : Group<Operation, Element>, 
    PartiallyInvertibleCommutativeMonoid<Operation, Element>,
    algebra::AbelianGroup<Operation, Element>
{};


// ========================
// Additive scalar concepts
// ========================


concept AdditiveMagma<typename Element>
  : Magma< math::add<Element>, Element >
{
    typename plus_assign_result_type;  
    plus_assign_result_type operator+=(Element& x, Element y);
    // requires std::Convertible<plus_assign_result_type, Element>;

    // Operator + is by default defined with +=
    typename addition_result_type;  
    addition_result_type operator+(Element x, Element y);
#if 0
    {
	Element tmp(x);
	return tmp += y;                      defaults NYS
    }
#endif 
    requires std::Convertible<addition_result_type, Element>;

    // Type consistency with Magma
    requires std::Convertible< addition_result_type,   
                            Magma< math::add<Element>, Element >::result_type>;

    // SameType requires more rigorous specializations on pure algebraic functors
    // requires std::SameType< addition_result_type, 
    // 	                    Magma< math::add<Element>, Element >::result_type>;

    axiom Consistency(math::add<Element> op, Element x, Element y)
    {
	op(x, y) == x + y;     
               
	// Consistency definition between + and += might change later
        x + y == x += y;
	// Element tmp = x; tmp+= y; tmp == x + y; not proposal-compliant
    }   
}


auto concept AdditiveCommutativeMagma<typename Element>
  : AdditiveMagma<Element>,
    CommutativeMagma< math::add<Element>, Element >
{};


auto concept AdditiveSemiGroup<typename Element>
  : AdditiveMagma<Element>, 
    SemiGroup< math::add<Element>, Element >
{};


// We really need only one of the additive concepts for the requirements, 
// the requirements of the other would be implied.
// Vice versa, to derive concept maps of nested concepts from
// concept maps of refined concepts, they are needed all.
auto concept AdditiveCommutativeSemiGroup<typename Element>
  : AdditiveSemiGroup<Element>,
    AdditiveCommutativeMagma<Element>,
    CommutativeSemiGroup< math::add<Element>, Element >
{};


concept AdditiveMonoid<typename Element>
  : AdditiveSemiGroup<Element>,
    Monoid< math::add<Element>, Element >
{
    Element zero(Element v);

    axiom Consistency (math::add<Element> op, Element x)
    {
	zero(x) == identity(op, x);
    }
};


// We really need only one of the additive concepts for the requirements, 
// the requirements of the other would be implied.
// Vice versa, to derive concept maps of nested concepts from
// concept maps of refined concepts, they are needed all.
auto concept AdditiveCommutativeMonoid<typename Element>
  : AdditiveMonoid<Element>,
    AdditiveCommutativeSemiGroup<Element>,
    CommutativeMonoid< math::add<Element>, Element >
{};


concept AdditivePartiallyInvertibleMonoid<typename Element>
  : AdditiveMonoid<Element>,
    PartiallyInvertibleMonoid< math::add<Element>, Element >
{
    typename minus_assign_result_type;  
    minus_assign_result_type operator-=(Element& x, Element y);
    // requires std::Convertible<minus_assign_result_type, Element>;
     
    // Operator - by default defined with -=
    typename subtraction_result_type;  
    subtraction_result_type operator-(Element& x, Element y);
#if 0
    {
	Element tmp(x);
	return tmp -= y;                      defaults NYS
    }
#endif 
    requires std::Convertible<subtraction_result_type, Element>;


    typename unary_result_type;  
    unary_result_type operator-(Element x);
#if 0
    {
	return zero(x) - x;      defaults NYS
    }
#endif 
    requires std::Convertible<unary_result_type, Element>;
    
    axiom Consistency(math::add<Element> op, Element x, Element y)
    {
	// consistency between additive and pure algebraic concept
	if ( is_invertible(op, y) )
	    op(x, inverse(op, y)) == x - y;            
	if ( is_invertible(op, y) )
	    inverse(op, y) == -y;                      

	// consistency between unary and binary -
	if ( is_invertible(op, x) )
	    identity(op, x) - x == -x;                 

	// Might change later
	if ( is_invertible(op, y) )
	    x - y == x -= y;                                                       
	// Element tmp = x; tmp-= y; tmp == x - y; not proposal-compliant
    }  

};


auto concept AdditivePartiallyInvertibleCommutativeMonoid<typename Element>
  : AdditivePartiallyInvertibleMonoid<Element>,
    AdditiveCommutativeMonoid<Element>, 
    PartiallyInvertibleCommutativeMonoid< math::add<Element>, Element >
{};



auto concept AdditiveGroup<typename Element>
  : AdditivePartiallyInvertibleMonoid<Element>,
    Group< math::add<Element>, Element >
{};


auto concept AdditiveAbelianGroup<typename Element>
  : AdditiveGroup<Element>,
    AdditiveCommutativeMonoid<Element>,
    AbelianGroup< math::add<Element>, Element >
{};


// ============================
// Multiplitive scalar concepts
// ============================


concept MultiplicativeMagma<typename Element>
  : Magma< math::mult<Element>, Element >
{
    typename mult_assign_result_type;  
    mult_assign_result_type operator*=(Element& x, Element y);
    // requires std::Convertible<mult_assign_result_type, Element>;

    // Operator * is by default defined with *=
    typename mult_result_type;  
    mult_result_type operator*(Element x, Element y);
#if 0
    {
	Element tmp(x);
	return tmp *= y;                      defaults NYS
    }
#endif 
    requires std::Convertible<mult_result_type, Element>;
    
    // Type consistency with Magma
    requires std::Convertible< mult_result_type,   
                            Magma< math::mult<Element>, Element >::result_type>;

    // SameType requires more rigorous specializations on pure algebraic functors
    // requires std::SameType< mult_result_type, 
    // 	                    Magma< math::mult<Element>, Element >::result_type>;


    axiom Consistency(math::mult<Element> op, Element x, Element y)
    {
	op(x, y) == x * y;                 
   
	// Consistency definition between * and *= might change later
        x * y == x *= y;
	// Element tmp = x; tmp*= y; tmp == x * y; not proposal-compliant
    }  

}


auto concept MultiplicativeSemiGroup<typename Element>
  : MultiplicativeMagma<Element>,
    SemiGroup< math::mult<Element>, Element >
{};


auto concept MultiplicativeCommutativeSemiGroup<typename Element>
  : MultiplicativeSemiGroup<Element>,
    CommutativeSemiGroup< math::mult<Element>, Element >
{};


concept MultiplicativeMonoid<typename Element>
  : MultiplicativeSemiGroup<Element>,
    Monoid< math::mult<Element>, Element >
{
    Element one(Element v);

    axiom Consistency (math::mult<Element> op, Element x)
    {
	one(x) == identity(op, x);
    }
};


auto concept MultiplicativeCommutativeMonoid<typename Element>
  : MultiplicativeMonoid<Element>,
    MultiplicativeCommutativeSemiGroup<Element>,
    CommutativeMonoid< math::mult<Element>, Element >
{};


concept MultiplicativePartiallyInvertibleMonoid<typename Element>
  : MultiplicativeMonoid<Element>,
    PartiallyInvertibleMonoid< math::mult<Element>, Element >
{
    typename divide_assign_result_type;  
    divide_assign_result_type operator/=(Element& x, Element y);
    // requires std::Convertible<divide_assign_result_type, Element>;
     
    // Operator / by default defined with /=
    typename division_result_type = Element;  
    division_result_type operator/(Element x, Element y);
#if 0
    {
	Element tmp(x);
	return tmp /= y;                      defaults NYS
    }
#endif 
    requires std::Convertible<division_result_type, Element>;
    
    axiom Consistency(math::mult<Element> op, Element x, Element y)
    {
	// consistency between multiplicative and pure algebraic concept
	if ( is_invertible(op, y) )
	    op(x, inverse(op, y)) == x / y;            

	// Consistency between / and /=, might change later
	if ( is_invertible(op, y) )
	    x / y == x /= y;              
	// Element tmp = x; tmp/= y; tmp == x / y; not proposal-compliant 
    }  
};
 

auto concept MultiplicativePartiallyInvertibleCommutativeMonoid<typename Element>
  : MultiplicativePartiallyInvertibleMonoid<Element>,
    MultiplicativeCommutativeMonoid<Element>,
    PartiallyInvertibleCommutativeMonoid< math::mult<Element>, Element >
{};
 

auto concept MultiplicativeGroup<typename Element>
  : MultiplicativeMonoid<Element>,
    Group< math::mult<Element>, Element >
{};


auto concept MultiplicativeAbelianGroup<typename Element>
  : MultiplicativeGroup<Element>,
    MultiplicativeCommutativeMonoid<Element>,
    AbelianGroup< math::mult<Element>, Element >
{};


// ======================================
// Algebraic concepts with two connectors
// ======================================

// -----------------
// Based on functors
// -----------------

// More generic, less handy to use

auto concept GenericRing<typename AddOp, typename MultOp, typename Element>
  : AbelianGroup<AddOp, Element>,
    SemiGroup<MultOp, Element>,
    algebra::Ring<AddOp, MultOp, Element>
{};


auto concept GenericCommutativeRing<typename AddOp, typename MultOp, typename Element>
  : GenericRing<AddOp, MultOp, Element>,
    CommutativeSemiGroup<MultOp, Element>
{};


auto concept GenericRingWithIdentity<typename AddOp, typename MultOp, typename Element>
  : GenericRing<AddOp, MultOp, Element>,
    Monoid<MultOp, Element>,
    algebra::RingWithIdentity<AddOp, MultOp, Element>
{};


auto concept GenericCommutativeRingWithIdentity<typename AddOp, typename MultOp, typename Element>
  : GenericRingWithIdentity<AddOp, MultOp, Element>,
    GenericCommutativeRing<AddOp, MultOp, Element>,
    CommutativeMonoid<MultOp, Element>
{};


// auto
concept GenericDivisionRing<typename AddOp, typename MultOp, typename Element>
  : GenericRingWithIdentity<AddOp, MultOp, Element>,
    algebra::DivisionRing<AddOp, MultOp, Element>
{
    requires std::Convertible<inverse_result_type, Element>;
};    


auto concept GenericField<typename AddOp, typename MultOp, typename Element>
  : GenericDivisionRing<AddOp, MultOp, Element>,
    GenericCommutativeRingWithIdentity<AddOp, MultOp, Element>,
    algebra::Field<AddOp, MultOp, Element>
{};


// ------------------
// Based on operators 
// ------------------

// Handier, less generic

// Alternative definitions use MultiplicativeMonoid<Element> for Ring
// and call such concepts Pseudo-Ring


auto concept Ring<typename Element>
  : AdditiveAbelianGroup<Element>,
    MultiplicativeSemiGroup<Element>,
    GenericRing<math::add<Element>, math::mult<Element>, Element>
{};


auto concept CommutativeRing<typename Element>
  : Ring<Element>,
    MultiplicativeCommutativeSemiGroup<Element>,
    GenericCommutativeRing<math::add<Element>, math::mult<Element>, Element>    
{};


auto concept RingWithIdentity<typename Element>
  : Ring<Element>,
    MultiplicativeMonoid<Element>,
    GenericRingWithIdentity<math::add<Element>, math::mult<Element>, Element>
{};
 

auto concept CommutativeRingWithIdentity<typename Element>
  : RingWithIdentity<Element>,
    CommutativeRing<Element>,
    MultiplicativeCommutativeMonoid<Element>,
    GenericCommutativeRingWithIdentity<math::add<Element>, math::mult<Element>, Element>
{};


concept DivisionRing<typename Element>
  : RingWithIdentity<Element>,
    MultiplicativePartiallyInvertibleMonoid<Element>, 
    GenericDivisionRing<math::add<Element>, math::mult<Element>, Element>
{
    axiom NonZeroDivisibility(Element x)
    {
	if (x != zero(x)) 
	    x / x == one(x);
    }
};    


auto concept Field<typename Element>
  : DivisionRing<Element>,
    CommutativeRingWithIdentity<Element>,
    GenericField<math::add<Element>, math::mult<Element>, Element>
{};


// ======================
// Miscellaneous concepts
// ======================

// that shall find a better place later


// EqualityComparable will have the != when defaults are supported
// At this point the following won't needed anymore
auto concept FullEqualityComparable<typename T, typename U = T>
{
  //requires std::EqualityComparable<T, U>;

    bool operator==(const T&, const U&);
    bool operator!=(const T&, const U&);
};

// Closure of EqualityComparable under a binary operation:
// That is, the result of this binary operation is also EqualityComparable
// with itself and with the operand type.
auto concept Closed2EqualityComparable<typename Operation, typename Element>
  : BinaryIsoFunction<Operation, Element>
{
    requires FullEqualityComparable<Element>;
    requires FullEqualityComparable< BinaryIsoFunction<Operation, Element>::result_type >;
    requires FullEqualityComparable< Element, BinaryIsoFunction<Operation, Element>::result_type >;
    requires FullEqualityComparable< BinaryIsoFunction<Operation, Element>::result_type, Element >;
};


// LessThanComparable will have the other operators when defaults are supported
// At this point the following won't needed anymore
auto concept FullLessThanComparable<typename T, typename U = T>
{
    bool operator<(const T&, const U&);
    bool operator<=(const T&, const U&);
    bool operator>(const T&, const U&);
    bool operator>=(const T&, const U&);
};


// Same for LessThanComparable
auto concept Closed2LessThanComparable<typename Operation, typename Element>
  : BinaryIsoFunction<Operation, Element>
{
    requires FullLessThanComparable<Element>;
    requires FullLessThanComparable< BinaryIsoFunction<Operation, Element>::result_type >;
    requires FullLessThanComparable< Element, BinaryIsoFunction<Operation, Element>::result_type >;
    requires FullLessThanComparable< BinaryIsoFunction<Operation, Element>::result_type, Element >;
};

#if 0
auto concept NumericOperatorResultConvertible<typename T>
  : AddableWithAssign<T>,
    SubtractableWithAssign<T>,
    MultiplicableWithAssign<T>,
    DivisibleWithAssign<T>
{
    requires std::Convertible< AddableWithAssign<T>::result_type, T>;
    requires std::Convertible< SubtractableWithAssign<T>::result_type, T>;
    requires std::Convertible< MultiplicableWithAssign<T>::result_type, T>;
    requires std::Convertible< DivisibleWithAssign<T>::result_type, T>;
}
#endif

auto concept AdditionResultConvertible<typename T>
{
    typename result_type;
    result_type operator+(T t, T u);
    requires std::Convertible<result_type, T>;

    typename result_type;
    result_type operator+=(T& t, T u);
    requires std::Convertible<result_type, T>;
};    


auto concept SubtractionResultConvertible<typename T>
{
    typename result_type;
    result_type operator-(T t, T u);
    requires std::Convertible<result_type, T>;

    typename result_type;
    result_type operator-=(T& t, T u);
    requires std::Convertible<result_type, T>;
};    

auto concept NumericOperatorResultConvertible<typename T>
  : AdditionResultConvertible<T>,
    SubtractionResultConvertible<T>
{};

// ====================
// Default Concept Maps
// ====================

#ifndef LA_NO_CONCEPT_MAPS

// ==============
// Integral Types
// ==============

template <typename T>
  requires IntrinsicSignedIntegral<T>
concept_map CommutativeRingWithIdentity<T> {}


template <typename T>
  requires IntrinsicUnsignedIntegral<T>
concept_map AdditiveCommutativeMonoid<T> {}

template <typename T>
  requires IntrinsicUnsignedIntegral<T>
concept_map MultiplicativeCommutativeMonoid<T> {}


// ====================
// Floating Point Types
// ====================

template <typename T>
  requires IntrinsicFloatingPoint<T>
concept_map Field<T> {}


template <typename T>
  requires IntrinsicFloatingPoint<T>
concept_map Field< std::complex<T> > {}


// ===========
// Min and Max
// ===========

// Draft version: defined generously unless there will be problems with some types

template <typename Element>
concept_map CommutativeMonoid< max<Element>, Element > 
{
    // Why do we need this?
    typedef Element identity_result_type;
}

template <typename Element>
concept_map CommutativeMonoid< min<Element>, Element >
{
    // Why do we need this?
    typedef Element identity_result_type;
}

#endif // LA_NO_CONCEPT_MAPS

// Here come some mathematical concepts to be defined later

/// Specify the semantic Behavior of natural numbers (TBD)
/** Mathematic properties are in most cases only approximated
    but not held exactly. Rigorous definition would impede
    usage of real processors. **/
concept NaturalNumber<typename T> {}

/// Specify the semantic Behavior of integral numbers (TBD)
/** Mathematic properties are in most cases only approximated
    but not held exactly. Rigorous definition would impede
    usage of real processors.
    It is arguable if this is really a refinement **/
concept IntegralNumber<typename T> : NaturalNumber<T> {}
    
/// Specify the semantic Behavior of complex numbers (TBD)
/** Mathematic properties are in most cases only approximated
    but not held exactly. Rigorous definition would impede
    usage of real processors. **/
concept ComplexNumber<typename T> {}

/// Specify the semantic Behavior of real numbers (TBD)
/** Mathematic properties are in most cases only approximated
    but not held exactly. Rigorous definition would impede
    usage of real processors. **/
concept RealNumber<typename T> : ComplexNumber<T> {}

#endif // __GXX_CONCEPTS__




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
    requires FullEqualityComparable<type>;
    requires FullLessThanComparable<type>;

    requires Field<type>;

    type sqrt(type);
    // typename sqrt_result;
    // sqrt_result sqrt(type);
    // requires std::Convertible<sqrt_result, type>;

    // using std::abs;
    type abs(T);
}

#else  // now without concepts

template <typename T>
struct Magnitude
{
    typename type = T;
};

template <typename T>
struct Magnitude<std::complex<T> >
{
    typedef T type;
}

template <typename T> struct RealMagnitude
  : public Magnitude<T>
{}

#endif  // __GXX_CONCEPTS__

// Type trait version both available with and w/o concepts (TBD: Macro finally :-( )
// For the moment everything is its own magnitude type, unless stated otherwise
template <typename T>
struct magnitude_type_trait
{
    typedef T type;
};

template <typename T>
struct magnitude_type_trait< std::complex<T> >
{
    typedef T type;
};


// =========================================
// Concepts for convenience (many from Rolf)
// =========================================


#ifdef __GXX_CONCEPTS__

//The following concepts Addable, Subtractable etc. differ from std::Addable, std::Subtractable 
//etc. in so far that no default for result_type is provided, thus allowing automated return type deduction

auto concept Addable<typename T, typename U = T>
{
    typename result_type;
    result_type operator+(const T& t, const U& u);
};
  
 
// Usually + and += are both defined
// + can be efficiently derived from += but not vice versa
auto concept AddableWithAssign<typename T, typename U = T>
{
    typename assign_result_type;  
    assign_result_type operator+=(T& x, U y);

    // Operator + is by default defined with +=
    typename result_type;  
    result_type operator+(T x, U y);
#if 0
    {
	// Default requires std::CopyConstructible, without default not needed
	Element tmp(x);                       
	return tmp += y;                      defaults NYS
    }
#endif 
};


auto concept Subtractable<typename T, typename U = T>
{
    typename result_type;
    result_type operator-(const T& t, const U& u);
};
  

// Usually - and -= are both defined
// - can be efficiently derived from -= but not vice versa
auto concept SubtractableWithAssign<typename T, typename U = T>
{
    typename assign_result_type;  
    assign_result_type operator-=(T& x, U y);

    // Operator - is by default defined with -=
    typename result_type;  
    result_type operator-(T x, U y);
#if 0
    {
	// Default requires std::CopyConstructible, without default not needed
	Element tmp(x);                       
	return tmp -= y;                      defaults NYS
    }
#endif 
};


auto concept Multiplicable<typename T, typename U = T>
{
    typename result_type;
    result_type operator*(const T& t, const U& u);
};


// Usually * and *= are both defined
// * can be efficiently derived from *= but not vice versa
auto concept MultiplicableWithAssign<typename T, typename U = T>
{
    typename assign_result_type;  
    assign_result_type operator*=(T& x, U y);

    // Operator * is by default defined with *=
    typename result_type;  
    result_type operator*(T x, U y);
#if 0
    {
	// Default requires std::CopyConstructible, without default not needed
	Element tmp(x);                       
	return tmp *= y;                      defaults NYS
    }
#endif 
};


auto concept Divisible<typename T, typename U = T>
{
    typename result_type;
    result_type operator / (const T&, const U&);
};


// Usually * and *= are both defined
// * can be efficiently derived from *= but not vice versa
auto concept DivisibleWithAssign<typename T, typename U = T>
{
    typename assign_result_type;  
    assign_result_type operator*=(T& x, U y);

    // Operator * is by default defined with *=
    typename result_type;  
    result_type operator*(T x, U y);
#if 0
    {
	// Default requires std::CopyConstructible, without default not needed
	Element tmp(x);                       
	return tmp *= y;                      defaults NYS
    }
#endif 
};


auto concept Transposable<typename T>
{
    typename result_type;
    result_type trans(T&);
};  


// Unary Negation -> Any suggestions for better names?! Is there a word as "negatable"?!
auto concept Negatable<typename S>
{
    typename result_type = S;
    result_type operator-(const S&);
};

// Or HasAbs?
using std::abs;
auto concept AbsApplicable<typename S>
{
    // There are better ways to define abs than the way it is done in std
    // Likely we replace the using one day
    typename result_type;
    result_type abs(const S&);
};


using std::conj;
auto concept HasConjugate<typename S>
{
    typename result_type;
    result_type conj(const S&);
};
  
  
// We need the following; might be placed somewhere else later
template <Float T>
concept_map HasConjugate<T> 
{ 
    typedef T result_type;
    result_type conj(const T& s) {return s;}
}



// Dot product to be defined:
auto concept Dottable<typename T, typename U = T>
{
    typename result_type = T;
    result_type dot(const T&t, const U& u);
};
    

auto concept OneNormApplicable<typename V> 
{
    typename result_type;
    result_type one_norm(const V&);
};


auto concept TwoNormApplicable<typename V> 
{
    typename result_type;
    result_type two_norm(const V&);
};


auto concept InfinityNormApplicable<typename V> 
{
    typename result_type;
    result_type inf_norm(const V&);
};




#endif  // __GXX_CONCEPTS__


} // namespace math



#endif // LA_CONCEPTS_INCLUDE
