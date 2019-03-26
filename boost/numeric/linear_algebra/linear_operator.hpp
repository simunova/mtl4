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

#ifndef MATH_LINEAR_OPERATOR_INCLUDE
#define MATH_LINEAR_OPERATOR_INCLUDE

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif

namespace math {

/** @addtogroup Concepts
 *  @{
 */

#ifdef __GXX_CONCEPTS__
    concept LinearOperator<typename Operator, typename VectorDomain, typename VectorImage>
    {
	requires VectorSpace<VectorDomain>;
	requires VectorSpace<VectorImage>;

	typename result_type;
	result_type operator* (Operator, VectorDomain);
	
	typename assign_type;
	assign_type operator= (VectorImage, result_type);

	// The following two requirements are subject to discussion
	typename plus_assign_type;
	plus_assign_type operator+= (VectorImage, result_type);
	
	typename minus_assign_type;
	minus_assign_type operator-= (VectorImage, result_type);

	axiom Addability(Operator A, VectorDomain x, VectorDomain y)
	{
	    A * (x + y) == A*x  + A*y;
	}

	// The two vector spaces must be scalable with the same scalar types
	axiom Scalability(Operator A, VectorSpace<VectorDomain>::scalar_type alpha, VectorDomain x)
	{
	    A * (alpha * x) == alpha * (A * x);
	}
    };
#else
    //! Concept LinearOperator
    /*!
        Linear operator from one vector space into another one.

        \param Operator The type of the operator, e.g., some matrix type
	\param VectorDomain The the type of a vector in the domain vector space
	\param VectorImage The the type of a vector in the image vector space
	
	\par Associated Types:
	- result_type
	- assign_type
	- plus_assign_type
	- minus_assign_type

	\par Requires:
	- VectorSpace < VectorDomain >
	- VectorSpace < VectorImage >

        \par Notation:
        <table summary="notation">
          <tr>
            <td>A</td>
	    <td>Object of type Operation</td>
	  </tr>
          <tr>
            <td>x, y</td>
	    <td>Objects of type VectorDomain</td>
          </tr>
          <tr>
            <td>u</td>
	    <td>Object of type VectorImage</td>
          </tr>
        </table>

        \par Valid Expressions:
        <table>
          <tr>
            <td>Assign product:</td>
	    <td>u= A * x</td>
	  </tr>
          <tr>
            <td>Add product:</td>
	    <td>u+= A * x</td>
	  </tr>
          <tr>
            <td>Subtract product:</td>
	    <td>u-= A * x</td>
	  </tr>
        </table>

        \invariant
        <table summary="invariants">
          <tr>
            <td>Addability</td>
	    <td>A * (x + y) == A*x + A*y</td>
          </tr>
          <tr>
            <td>Scalability</td>
	    <td>alpha * (A * x) == A * (alpha * x)</td>
          </tr>
        </table>
	
	\note
	-# Using matrix vector products in arbitrary expressions requires
	   storing it in temporary objects to avoid redundant computation.
	   On the other hand, it is not always obvious to choose an appropriate
	   type for such temporary depending on arbitrary operator and vector types.
	   Using the products directly in assignments allows implementation without
	   temporaries, e.g., by calling a function mult(A, x, u) internally.
     */
    template <typename Operator, typename VectorDomain, typename VectorImage>
    struct LinearOperator
    {
	/// Associated type: result of multiplication; automatically deducted
	typedef associated_type result_type;
	/// Multiplication of linear operator with vector
	result_type operator* (Operator, VectorDomain);

	/// Associated type: return type of assigning product to vector.
	/** Automatically deducted. Using expression templates it can be different from VectorImage& */
	typedef associated_type assign_type;
	/// Product must be assignable
	assign_type operator= (VectorImage, result_type);

	// The following two requirements are subject to discussion
	/// Associated type: return type of incrementally assigning product to vector.
	/** Automatically deducted. Using expression templates it can be different from VectorImage& */
	typedef associated_type plus_assign_type;
	/// Product must be assignable with increment
	plus_assign_type operator+= (VectorImage, result_type);
	
	// The following two requirements are subject to discussion
	/// Associated type: return type of decrementally assigning product to vector.
	/** Automatically deducted. Using expression templates it can be different from VectorImage& */
	typedef associated_type minus_assign_type;
	/// Product must be assignable with decrement
	minus_assign_type operator+= (VectorImage, result_type);

	/// Invariant: the linear projection of a sum is the sum of the linear projections
	axiom Addability(Operator A, VectorDomain x, VectorDomain y)
	{
	    A * (x + y) == A*x  + A*y;
	}

	/// Invariant: the linear projection of a scaled vector is the scaling of the vector's linear projections
	axiom Scalability(Operator A, VectorSpace<VectorDomain>::scalar_type alpha, VectorDomain x)
	{
	    A * (alpha * x) == alpha * (A * x);
	}	
    };
#endif


#ifdef __GXX_CONCEPTS__
    concept SelfAdjointOperator<typename Operator, typename VectorDomain, typename VectorImage>
      : LinearOperator<Operator, VectorDomain, VectorImage>
    {};
#else
    //! Concept SelfAdjointOperator
    /*!
        

        \param Operator The type of the operator, e.g., some matrix type
	\param VectorDomain The the type of a vector in the domain vector space
	\param VectorImage The the type of a vector in the image vector space
	
        \par Refinement of:
	- LinearOperator <Operator, VectorDomain, VectorImage>
    */
    template <typename Operator, typename VectorDomain, typename VectorImage>
    struct SelfAdjointOperator
      : LinearOperator<Operator, VectorDomain, VectorImage>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept RealOperator<typename Operator, typename VectorDomain, typename VectorImage>
      : LinearOperator<Operator, VectorDomain, VectorImage>
    {};
#else
    //! Concept RealOperator
    /*!
        

        \param Operator The type of the operator, e.g., some matrix type
	\param VectorDomain The the type of a vector in the domain vector space
	\param VectorImage The the type of a vector in the image vector space
	
        \par Refinement of:
	- LinearOperator <Operator, VectorDomain, VectorImage>
    */
    template <typename Operator, typename VectorDomain, typename VectorImage>
    struct RealOperator
      : LinearOperator<Operator, VectorDomain, VectorImage>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept SymmetricOperator<typename Operator, typename VectorDomain, typename VectorImage>
      : SelfAdjointOperator<Operator, VectorDomain, VectorImage>,
        RealOperator<Operator, VectorDomain, VectorImage>
    {};
#else
    //! Concept SymmetricOperator
    /*!
        

        \param Operator The type of the operator, e.g., some matrix type
	\param VectorDomain The the type of a vector in the domain vector space
	\param VectorImage The the type of a vector in the image vector space
	
        \par Refinement of:
	- SelfAdjointOperator <Operator, VectorDomain, VectorImage>
	- RealOperator <Operator, VectorDomain, VectorImage>
    */
    template <typename Operator, typename VectorDomain, typename VectorImage>
    struct SymmetricOperator
      : SelfAdjointOperator<Operator, VectorDomain, VectorImage>,
        RealOperator<Operator, VectorDomain, VectorImage>
    {};
#endif

/*@}*/ // end of group Concepts

} // namespace math

#endif // MATH_LINEAR_OPERATOR_INCLUDE
