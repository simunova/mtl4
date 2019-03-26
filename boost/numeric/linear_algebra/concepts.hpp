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

#ifndef MTL_NEW_ALGEBRAIC_CONCEPTS_INCLUDE
#define MTL_NEW_ALGEBRAIC_CONCEPTS_INCLUDE

#include <concepts>
#include <boost/numeric/linear_algebra/intrinsic_concept_maps.hpp>
#include <boost/numeric/linear_algebra/operators.hpp>


namespace math {

concept Commutative<typename Operation, typename Element>
  : std::Callable2<Op, Element, Element>
{
    axiom Commutativity(Operation op, Element x, Element y)
    {
	op(x, y) == op(y, x); 
    }   
};


concept SemiGroup<typename Operation, typename Element>
  : std::Callable2<Op, Element, Element>
{
    axiom Associativity(Operation op, Element x, Element y, Element z)
    {
	op(x, op(y, z)) == op(op(x, y), z); 
    }
};


concept Monoid<typename Operation, typename Element>
  : SemiGroup<Operation, Element> 
{
    typename identity_result_type;
    identity_result_type identity(Operation, Element);
    
    axiom Neutrality(Operation op, Element x)
    {
	op( x, identity(op, x) ) == x;
	op( identity(op, x), x ) == x;
    }
};


auto concept Inversion<typename Operation, typename Element>
{
    typename result_type;
    result_type inverse(Operation, Element);
    
};


concept PIMonoid<typename Operation, typename Element>
  : Monoid<Operation, Element>, 
    Inversion<Operation, Element>
{
    bool is_invertible(Operation, Element);
    
    requires std::Convertible<Inversion<Operation, Element>::result_type, Element>;

    axiom Invertibility(Operation op, Element x)
    {
	// Only for invertible elements:
	if (is_invertible(op, x))
	    op( x, inverse(op, x) ) == identity(op, x); 
	if ( is_invertible(op, x) )
	    op( inverse(op, x), x ) == identity(op, x); 
    }
}

#if 0 
    // Alternative approach to convert the result of inversion to Element
    // Unfortunately, this doesn't compile
    template <typename Operation, typename Element>
        requires PIMonoid<Operation, Element>
    concept_map PIMonoid<Operation, Inversion<Operation, Element>::result_type> {}
#endif


concept Group<typename Operation, typename Element>
  : PIMonoid<Operation, Element>
{
    bool is_invertible(Operation, Element) { return true; }
    
    // Just in case somebody redefines is_invertible
    axiom AlwaysInvertible(Operation op, Element x)
    {
	is_invertible(op, x);
    }

    // In fact this is implied by AlwaysInvertible and inherited Invertibility axiom
    // Maybe remove
    axiom GlobalInvertibility(Operation op, Element x)
    {
	op( x, inverse(op, x) ) == identity(op, x);
	op( inverse(op, x), x ) == identity(op, x);
    }
};


auto concept AbelianGroup<typename Operation, typename Element>
  : Group<Operation, Element>, Commutative<Operation, Element>
{};
    

// =======================
// Operator-based concepts
// =======================


concept Additive<typename Element>
  : std::HasPlus<Element>
{
    typename plus_assign_result_type;  
    plus_assign_result_type operator+=(Element& x, Element y)
    {
	x= x + y; return x;
    }
    
    requires std::Convertible<plus_assign_result_type, Element&>;
    
    // Do we need the opposite conversion too?
    // This line produces a compiler error
    // requires std::Convertible<add<Element>::result_type,
    //                           std::HasPlus<Element>::result_type>;

    axiom Consistency(add<Element> op, Element x, Element y)
    {
	op(x, y) == x + y;
	op(x, y) == (x += y, x);
    }
}

    
auto concept AdditiveCommutative<typename Element>
  : Additive<Element>,
    Commutative< add<Element>, Element >
{}

    
auto concept AdditiveSemiGroup<typename Element>
  : Additive<Element>,
    SemiGroup< add<Element>, Element >
{}

    
concept AdditiveMonoid<typename Element>
  : AdditiveSemiGroup<Element>,
    Monoid< add<Element>, Element >
{
    Element zero(Element x)
    {
	return identity(add<Element>(), x);
    }
    
    // If we don't use the default definition
    axiom IdentityConsistency (add<Element> op, Element x)
    {
	zero(x) == identity(op, x);
    }
};

#ifndef CONCEPTS_WITHOUT_OVERLOADED_REQUIREMENTS  // Uncompilable due to error in compiler
concept AdditivePIMonoid<typename Element>
  : std::HasMinus<Element>, AdditiveMonoid<Element>, 
    PIMonoid< add<Element>, Element >
{
    typename minus_assign_result_type;  
    minus_assign_result_type operator-=(Element& x, Element y)
    {
	x= x - y; return x;
    }
    
    requires std::Convertible<minus_assign_result_type, Element&>;
    
    typename unary_result_type;  
    unary_result_type operator-(Element x)
    {
	return zero(x) - x;
    }
    
    axiom InverseConsistency(add<Element> op, Element x, Element y)
    {
	// consistency between additive and functor concept
	if ( is_invertible(op, y) )
	    op(x, inverse(op, y)) == x - y;
	if ( is_invertible(op, y) )
	    op(x, y) == (x -= y, x);
	
	// consistency of unary inversion
	if ( is_invertible(op, y) )
	    inverse(op, y) == -y;                      

	// consistency between unary and binary -
	if ( is_invertible(op, x) )
	    identity(op, x) - x == -x;                 
    }
}


auto concept AdditiveGroup<typename Element>
  : AdditivePIMonoid<Element>,
    Group< add<Element>, Element >
{};


auto concept AdditiveAbelianGroup<typename Element>
  : AdditiveGroup<Element>,
    Commutative< add<Element>, Element >
{}

#endif


concept Multiplicative<typename Element>
  : std::HasMultiply<Element>
{
    typename times_assign_result_type;  
    times_assign_result_type operator*=(Element& x, Element y)
    {
	x= x * y; return x;
    }
	
    requires std::Convertible<times_assign_result_type, Element&>;
    
    // Do we need the opposite conversion too?
    // This line produces a compiler error
    // requires std::Convertible<mult<Element>::result_type,
    //                           std::HasMultiply<Element>::result_type>;

    axiom Consistency(mult<Element> op, Element x, Element y)
    {
	op(x, y) == x * y;
	op(x, y) == (x *= y, x);
    }
}


auto concept MultiplicativeCommutative<typename Element>
  : Multiplicative<Element>,
    Commutative< mult<Element>, Element >
{}

    
auto concept MultiplicativeSemiGroup<typename Element>
  : Multiplicative<Element>,
    SemiGroup< mult<Element>, Element >
{}


concept MultiplicativeMonoid<typename Element>
  : MultiplicativeSemiGroup<Element>,
    Monoid< mult<Element>, Element >
{
    Element one(Element x)
    {
	return identity(mult<Element>(), x);
    }
	
    // If we don't use the default definition
    axiom IdentityConsistency (math::mult<Element> op, Element x)
    {
	one(x) == identity(op, x);
    }
};

#ifndef CONCEPTS_WITHOUT_OVERLOADED_REQUIREMENTS  // Uncompilable due to error in compiler
concept MultiplicativePIMonoid<typename Element>
  : std::HasDivide<Element>, MultiplicativeMonoid<Element>,
    PIMonoid< mult<Element>, Element >
{
    typename divide_assign_result_type;  
    divide_assign_result_type operator/=(Element& x, Element y)
    {
	x= x / y; return x;
    }
    
    requires std::Convertible<divide_assign_result_type, Element&>;
    
    axiom InverseConsistency(mult<Element> op, Element x, Element y)
    {
	// consistency between multiplicative and functor concept
	if ( is_invertible(op, y) )
	    op(x, inverse(op, y)) == x / y;
	if ( is_invertible(op, y) )
	    op(x, y) == (x /= y, x);
    }
}


auto concept MultiplicativeGroup<typename Element>
  : MultiplicativePIMonoid<Element>,
    Group< mult<Element>, Element >
{};


auto concept MultiplicativeAbelianGroup<typename Element>
  : MultiplicativeGroup<Element>,
    Commutative< mult<Element>, Element >
{}


// ==========================
// Concepts with 2 operations
// ==========================



concept Distributive<typename AddOp, typename MultOp, typename Element>
{
    axiom Distributivity(AddOp add, MultOp mult, Element x, Element y, Element z)
    {
	// From left
	mult(x, add(y, z)) == add(mult(x, y), mult(x, z));
	// from right
	mult(add(x, y), z) == add(mult(x, z), mult(y, z));
    }
}


auto concept Ring<typename AddOp, typename MultOp, typename Element>
  : AbelianGroup<AddOp, Element>,
    SemiGroup<MultOp, Element>,
    Distributive<AddOp, MultOp, Element>
{}


auto concept RingWithIdentity<typename AddOp, typename MultOp, typename Element>
  : Ring<AddOp, MultOp, Element>,
    Monoid<MultOp, Element>
{}


concept DivisionRing<typename AddOp, typename MultOp, typename Element>
  : RingWithIdentity<AddOp, MultOp, Element>,
    Inversion<MultOp, Element>
{
    // 0 != 1, otherwise trivial
    axiom ZeroIsDifferentFromOne(AddOp add, MultOp mult, Element x)
    {
	identity(add, x) != identity(mult, x);       
    }
    
    // Non-zero divisibility from left and from right
    axiom NonZeroDivisibility(AddOp add, MultOp mult, Element x)
    {
	if (x != identity(add, x))
	    mult(inverse(mult, x), x) == identity(mult, x);
	if (x != identity(add, x))
	    mult(x, inverse(mult, x)) == identity(mult, x);
    }
}    


auto concept Field<typename AddOp, typename MultOp, typename Element>
  : DivisionRing<AddOp, MultOp, Element>,
    Commutative<MultOp, Element>
{}


auto concept OperatorRing<typename Element>
  : AdditiveAbelianGroup<Element>,
    MultiplicativeSemiGroup<Element>,
    Ring<add<Element>, mult<Element>, Element>
{}

        
auto concept OperatorRingWithIdentity<typename Element>
  : OperatorRing<Element>,
    MultiplicativeMonoid<Element>,
    RingWithIdentity<add<Element>, mult<Element>, Element>
{}
       
 
auto concept OperatorDivisionRing<typename Element>
  : OperatorRingWithIdentity<Element>,
    MultiplicativePIMonoid<Element>, 
    DivisionRing<add<Element>, mult<Element>, Element>
{}    


auto concept OperatorField<typename Element>
  : OperatorDivisionRing<Element>,
    Field<add<Element>, mult<Element>, Element>
{}
        


#endif

concept IntrinsicType<typename T> {}

concept IntrinsicArithmetic<typename T> : IntrinsicType<T> {}

concept IntrinsicIntegral<typename T> : IntrinsicArithmetic<T> {}

concept IntrinsicSignedIntegral<typename T> 
  : std::SignedIntegralLike<T>,
    IntrinsicIntegral<T>
{}

concept IntrinsicUnsignedIntegral<typename T> 
  : std::UnsignedIntegralLike<T>,
    IntrinsicIntegral<T>
{}

concept IntrinsicFloatingPoint<typename T>
  : std::FloatingPointLike<T>,
    IntrinsicArithmetic<T>
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


#ifndef CONCEPTS_WITHOUT_OVERLOADED_REQUIREMENTS

// ====================
// Default Concept Maps
// ====================

// ==============
// Arithmetic
// ==============

// ----------------
// Signed integrals
// ----------------

template <typename T>
  requires IntrinsicSignedIntegral<T>
concept_map OperatorRingWithIdentity<T> {}

template <typename T>
  requires IntrinsicSignedIntegral<T>
concept_map MultiplicativeCommutative<T> {}

// ------------------
// Unsigned integrals
// ------------------


template <typename T>
  requires IntrinsicUnsignedIntegral<T>
concept_map AdditiveCommutative<T> {}

template <typename T>
  requires IntrinsicUnsignedIntegral<T>
concept_map AdditiveMonoid<T> {}

template <typename T>
  requires IntrinsicUnsignedIntegral<T>
concept_map MultiplicativeCommutative<T> {}

template <typename T>
  requires IntrinsicUnsignedIntegral<T>
concept_map MultiplicativeMonoid<T> {}

// ---------------
// Floationg Point
// ---------------


template <typename T>
  requires IntrinsicFloatingPoint<T>
concept_map Field<T> {}

template <typename T>
  requires IntrinsicFloatingPoint<T>
concept_map Field< std::complex<T> > {}


// ===========
// Min and Max
// ===========


template <typename T>
  requires IntrinsicArithmetic<T>
concept_map Commutative< max<T>, T > {}

template <typename T>
  requires IntrinsicArithmetic<T>
concept_map Monoid< max<T>, T > {}

template <typename T>
  requires IntrinsicArithmetic<T>
concept_map Commutative< min<T>, T > {}

template <typename T>
  requires IntrinsicArithmetic<T>
concept_map Monoid< min<T>, T > {}


// ==========
// And and Or
// ==========

template <typename T>
  requires Intrinsic<T> && std::HasLogicalAnd<T>
concept_map Commutative< std::logical_and<T>, T > {}

template <typename T>
  requires Intrinsic<T> && std::HasLogicalAnd<T>
concept_map Monoid< std::logical_and<T>, T > {}

template <typename T>
  requires Intrinsic<T> && std::HasLogicalOr<T>
concept_map Commutative< std::logical_or<T>, T > {}

template <typename T>
  requires Intrinsic<T> && std::HasLogicalOr<T>
concept_map Monoid< std::logical_or<T>, T > {}

template <typename T>
  requires Intrinsic<T> && std::HasLogicalAnd<T> && std::HasLogicalOr<T>
concept_map Distributive<std::logical_and<T>, std::logical_or<T>, T> {}

template <typename T>
  requires Intrinsic<T> && std::HasLogicalAnd<T> && std::HasLogicalOr<T>
concept_map Distributive<std::logical_or<T>, std::logical_and<T>, T> {}


// ==================
// Bitwise operations
// ==================

// not yet defined

template <typename T>
  requires IntrinsicIntegral<T>
concept_map Commutative< bit_and<T>, T > {}

template <typename T>
  requires IntrinsicIntegral<T>
concept_map Monoid< bit_and<T>, T > {}

template <typename T>
  requires IntrinsicIntegral<T>
concept_map Commutative< bit_or<T>, T > {}

template <typename T>
  requires IntrinsicIntegral<T>
concept_map Monoid< bit_or<T>, T > {}

template <typename T>
  requires IntrinsicIntegral<T>
concept_map Distributive<bit_and<T>, bit_or<T>, T> {}

template <typename T>
  requires IntrinsicIntegral<T> 
concept_map Distributive<bit_or<T>, bit_and<T>, T> {}

template <typename T>
  requires IntrinsicIntegral<T>
concept_map Commutative< bit_xor<T>, T > {}

template <typename T>
  requires IntrinsicIntegral<T>
concept_map SemiGroup< bit_xor<T>, T > {}

// ====================
// String concatenation
// ====================

concept_map AdditiveMonoid<std::string> {}


#endif



} // namespace math

#endif // MTL_NEW_ALGEBRAIC_CONCEPTS_INCLUDE
