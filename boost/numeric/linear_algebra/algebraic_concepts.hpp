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

#ifndef LA_ALGEBRAIC_CONCEPTS_DOC_INCLUDE
#define LA_ALGEBRAIC_CONCEPTS_DOC_INCLUDE

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>

/// Namespace for purely algebraic concepts
namespace algebra {

/** @addtogroup Concepts
 *  @{
 */

#ifndef __GXX_CONCEPTS__
    //! Concept Commutative
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined   

        \par Notation:
        <table summary="notation">
          <tr>
            <td>op</td>
    	<td>Object of type Operation</td>
          </tr>
          <tr>
            <td>x, y</td>
    	<td>Objects of type Element</td>
          </tr>
        </table>
        \invariant
        <table summary="invariants">
          <tr>
            <td>Commutativity</td>
    	<td>op(x, y) == op(y, x)</td>
          </tr>
        </table>
     */
    template <typename Operation, typename Element>
    struct Commutative 
    {};
#else 
    concept Commutative<typename Operation, typename Element>
    {
	axiom Commutativity(Operation op, Element x, Element y)
	{
	    op(x, y) == op(y, x); 
	}   
    };
#endif


#ifndef __GXX_CONCEPTS__
    //! Concept Associative
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined   

        \par Notation:
        <table summary="notation">
          <tr>
            <td>op</td>
    	<td>Object of type Operation</td>
          </tr>
          <tr>
            <td>x, y, z</td>
    	<td>Objects of type Element</td>
          </tr>
        </table>
        \invariant
        <table summary="invariants">
          <tr>
            <td>Associativity</td>
    	<td>op(x, op(y, z)) == op(op(x, y), z)</td>
          </tr>
        </table>
     */
    template <typename Operation, typename Element>
    struct Associative
    {};
#else
    concept Associative<typename Operation, typename Element>
    {
        axiom Associativity(Operation op, Element x, Element y, Element z)
        {
	    op(x, op(y, z)) == op(op(x, y), z); 
        }
    };
#endif


#ifndef __GXX_CONCEPTS__
    //! Concept SemiGroup
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined   

        \note
        -# The algebraic concept SemiGroup only requires associativity and is identical with the concept Associative.
     */
    template <typename Operation, typename Element>
    struct SemiGroup
        : Associative<Operation, Element>
    {};
#else
    auto concept SemiGroup<typename Operation, typename Element>
      : Associative<Operation, Element>
    {};
#endif


#ifndef __GXX_CONCEPTS__
    //! Concept Monoid
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- SemiGroup
        \par Notation:
        <table summary="notation">
          <tr>
            <td>op</td>
    	<td>Object of type Operation</td>
          </tr>
          <tr>
            <td>x</td>
    	<td>Object of type Element</td>
          </tr>
        </table>
        \invariant
        <table summary="invariants">
          <tr>
            <td>Neutrality from right</td>
	    <td>op( x, identity(op, x) ) == x</td>
          </tr>
          <tr>
            <td>Neutrality from left</td>
	    <td>op( identity(op, x), x ) == x</td>
          </tr>
        </table>
     */
    template <typename Operation, typename Element>
    struct Monoid
      : SemiGroup<Operation, Element> 
    {
	/// Associated type; if not defined in concept_map automatically detected as result of identity
        typedef associated_type identity_result_type; 
        identity_result_type identity(Operation, Element); ///< Identity element of Operation
    };
#else
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
#endif

#ifdef __GXX_CONCEPTS__
    auto concept Inversion<typename Operation, typename Element>
    {
        typename inverse_result_type;
        inverse_result_type inverse(Operation, Element);
     
    };
#else
    //! Concept Inversion
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined  

	\par Associated Types:
	- inverse_result_type
	\par Valid Expressions:
	- inverse(op, x);
     */
    template <typename Operation, typename Element>
    struct Inversion
    {
	/// Associated type; if not defined in concept_map automatically detected as result of inverse
        typedef associated_type inverse_result_type;

	/// Returns inverse of \p x regarding operation \p op
        inverse_result_type inverse(Operation op, Element x);
    };
#endif


#ifdef __GXX_CONCEPTS__
    concept Group<typename Operation, typename Element>
      : Monoid<Operation, Element>, Inversion<Operation, Element>
    {
        axiom Inversion(Operation op, Element x)
        {
	    op( x, inverse(op, x) ) == identity(op, x);
	    op( inverse(op, x), x ) == identity(op, x);
        }
    };
#else
    //! Concept Group
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- Monoid
	- Inversion
        \par Notation:
        <table summary="notation">
          <tr>
            <td>op</td>
    	<td>Object of type Operation</td>
          </tr>
          <tr>
            <td>x</td>
    	<td>Object of type Element</td>
          </tr>
        </table>
        \invariant
        <table summary="invariants">
          <tr>
            <td>Inverse from right</td>
    	<td>op( x, inverse(op, x) ) == identity(op, x)</td>
          </tr>
          <tr>
            <td>Inverse from left</td>
    	<td>op( inverse(op, x), x ) == identity(op, x)</td>
          </tr>
        </table>
     */
    template <typename Operation, typename Element>
    struct Group
      : Monoid<Operation, Element>,
	Inversion<Operation, Element>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    auto concept AbelianGroup<typename Operation, typename Element>
      : Group<Operation, Element>, Commutative<Operation, Element>
    {};
#else
    //! Concept AbelianGroup
    /*!
        \param Operation A functor implementing a binary operation
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- Group
	- Commutative
     */
    template <typename Operation, typename Element>
    struct AbelianGroup
      : Group<Operation, Element>,
	Commutative<Operation, Element>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept Distributive<typename AddOp, typename MultOp, typename Element>
    {
        axiom Distributivity(AddOp add, MultOp mult, Element x, Element y, Element z)
        {
	    // From left
	    mult(x, add(y, z)) == add(mult(x, y), mult(x, z));
	    // z right
	    mult(add(x, y), z) == add(mult(x, z), mult(y, z));
        }
    };
#else
    //! Concept Distributive
    /*!
        \param AddOp A functor implementing a binary operation representing addition
        \param MultOp A functor implementing a binary operation representing multiplication
        \param Element The type upon the binary operation is defined   

        \par Notation:
        <table summary="notation">
          <tr>
            <td>add</td>
	    <td>Object of type AddOp</td>
          </tr>
          <tr>
            <td>mult</td>
	    <td>Object of type Multop</td>
          </tr>
          <tr>
            <td>x, y, z</td>
	    <td>Objects of type Element</td>
          </tr>
        </table>
        \invariant
        <table summary="invariants">
          <tr>
            <td>Distributivity from left</td>
	    <td>mult(x, add(y, z)) == add(mult(x, y), mult(x, z))</td>
          </tr>
          <tr>
            <td>Distributivity from right</td>
	    <td>mult(add(x, y), z) == add(mult(x, z), mult(y, z))</td>
          </tr>
        </table>
     */    
    template <typename AddOp, typename MultOp, typename Element>
    struct Distributive
    {};
#endif


#ifdef __GXX_CONCEPTS__
    auto concept Ring<typename AddOp, typename MultOp, typename Element>
      : AbelianGroup<AddOp, Element>,
        SemiGroup<MultOp, Element>,
        Distributive<AddOp, MultOp, Element>
    {};
#else
    //! Concept Ring
    /*!
        \param AddOp A functor implementing a binary operation representing addition
        \param MultOp A functor implementing a binary operation representing multiplication
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- AbelianGroup <MultOp, Element>
	- SemiGroup <MultOp, Element>
        - Distributive <AddOp, MultOp, Element>
     */
    template <typename AddOp, typename MultOp, typename Element>
    struct Ring
      : AbelianGroup<AddOp, Element>,
        SemiGroup<MultOp, Element>,
        Distributive<AddOp, MultOp, Element>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    auto concept RingWithIdentity<typename AddOp, typename MultOp, typename Element>
      : Ring<AddOp, MultOp, Element>,
        Monoid<MultOp, Element>
    {};
#else
    //! Concept RingWithIdentity
    /*!
        \param AddOp A functor implementing a binary operation representing addition
        \param MultOp A functor implementing a binary operation representing multiplication
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- Ring <AddOp, MultOp, Element>
        - Monoid <MultOp, Element>
     */
    template <typename AddOp, typename MultOp, typename Element>
    struct RingWithIdentity
      : Ring<AddOp, MultOp, Element>,
        Monoid<MultOp, Element>
    {};
#endif


#ifdef __GXX_CONCEPTS__
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
    };    
#else
    //! Concept DivisionRing
    /*!
        \param AddOp A functor implementing a binary operation representing addition
        \param MultOp A functor implementing a binary operation representing multiplication
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- RingWithIdentity <AddOp, MultOp, Element>
        - Inversion <MultOp, Element>

        \par Notation:
        <table summary="notation">
          <tr>
            <td>add</td>
	    <td>Object of type AddOp</td>
          </tr>
          <tr>
            <td>mult</td>
	    <td>Object of type Multop</td>
          </tr>
          <tr>
            <td>x, y, z</td>
	    <td>Objects of type Element</td>
          </tr>
        </table>

        \invariant
        <table summary="invariants">
          <tr>
            <td>Non-zero divisibility from left</td>
	    <td>mult(inverse(mult, x), x) == identity(mult, x)</td>
	    <td>if x != identity(add, x)</td>
          </tr>
          <tr>
            <td>Non-zero divisibility from right</td>
	    <td>mult(x, inverse(mult, x)) == identity(mult, x)</td>
	    <td>if x != identity(add, x)</td>
          </tr>
          <tr>
            <td>Zero is different from one</td>
	    <td>identity(add, x) != identity(mult, x)</td>
	    <td></td>
          </tr>
	</table>

	\note
	-# Zero and one can be theoretically identical in a DivisionRing.  However,
	   this implies that there is only one element x in the Ring with x + x = x and 
	   x * x = x (which is actually even a Field). 
	   Because this structure has no practical value we exclude it from 
	   consideration.
     */
    template <typename AddOp, typename MultOp, typename Element>
    struct DivisionRing
      : RingWithIdentity<AddOp, MultOp, Element>,
        Inversion<MultOp, Element>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    // SkewField is defined as synonym for DivisionRing
    auto concept SkewField<typename AddOp, typename MultOp, typename Element>
      : DivisionRing<AddOp, MultOp, Element>
    {};
#else
    //! Concept SkewField
    /*!
        \param AddOp A functor implementing a binary operation representing addition
        \param MultOp A functor implementing a binary operation representing multiplication
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- DivisionRing <AddOp, MultOp, Element>
	\note
	- Because the refinement of DivisionRing to SkewField is automatic the two concepts
	  are identical.
     */
    template <typename AddOp, typename MultOp, typename Element>
    struct SkewField
      : DivisionRing<AddOp, MultOp, Element>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    auto concept Field<typename AddOp, typename MultOp, typename Element>
      : DivisionRing<AddOp, MultOp, Element>,
        Commutative<MultOp, Element>
    {};
#else
    //! Concept Field
    /*!
        \param AddOp A functor implementing a binary operation representing addition
        \param MultOp A functor implementing a binary operation representing multiplication
        \param Element The type upon the binary operation is defined   

        \par Refinement of:
	- DivisionRing <AddOp, MultOp, Element>
	- Commutative <MultOp, Element>
     */
    template <typename AddOp, typename MultOp, typename Element>
    struct Field
      : DivisionRing<AddOp, MultOp, Element>,
        Commutative<MultOp, Element>
    {};
#endif

/*@}*/ // end of group Concepts

} // algebra

#endif // LA_ALGEBRAIC_CONCEPTS_DOC_INCLUDE
