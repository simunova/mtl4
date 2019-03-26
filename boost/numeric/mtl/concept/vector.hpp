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

#ifndef MTL_VECTOR_CONCEPTS_INCLUDE
#define MTL_VECTOR_CONCEPTS_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif

namespace mtl {

/** @addtogroup Concepts
 *  @{
 */

#ifdef __GXX_CONCEPTS__
    concept Vector<typename T>
      : AlgebraicCollection<T>
    {
	const_reference T::operator() (size_type index) const;

	const_reference T::operator[] (size_type index) const;

	size_type nnz(T); // maybe into AlgebraicCollection?
    };
#else
    /// Concept Vector
    /**
        \par Refinement of:
	- AlgebraicCollection < T >
	\par Notation:
	- X is a type that models Vector
	- v is an object of type X
	- r are objects of size_type
	\par Valid expressions:
	- Element access: \n v(r) \n Return Type: const_reference 
	  \n Semantics: Element in row \p r and column \p c
	- Element access: \n v[r] \n Equivalent to v(r)
	\invariant
	- Either num_cols(v), in case of a column vector, or num_rows(v),
	  in case of a row vector, must be 1! Otherwise it would be
	  a matrix. 
	\par Models:
       	- dense_vector
	\note
	-# If it would become (extremely) important to support 1D C arrays as Vector, 
	  it might be necessary to drop the requirement
	  of element access by v(r).
     */ 
    template <typename T>
    struct Vector
	: public AlgebraicCollection<T>
    {
	/// Element access
	const_reference T::operator() (size_type index) const;

	/// Element access
	const_reference T::operator[] (size_type index) const;
    };
#endif


#ifdef __GXX_CONCEPTS__
    concept MutableVector<typename T>
      : Vector<T>, 
        MutableCollection<T>
    {
	reference T::operator() (size_type index);

	reference T::operator[] (size_type index);
    };
#else
    /// Concept MutableVector
    /**
        \par Refinement of:
	- Vector < T >
	- MutableCollection < T >
	\par Notation:
	- X is a type that models Vector
	- v is an object of type X
	- r are objects of size_type
	\par Valid expressions:
	- Element access: \n v(r) \n Return Type: reference 
	  \n Semantics: Element in row \p r for a column vector or in column \p c for a row vector
	- Element access: \n v[r] \n Equivalent to v(r)
	\par Models:
       	- dense_vector
	\note
	-# If it would become (extremely) important to support 1D C arrays as Vector, 
	  it might be necessary to drop the requirement
	  of element access by v(r).
     */ 
    template <typename T>
    struct MutableVector
      : public Vector<T>, 
        public MutableCollection<T>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept ConstantSizeVector<typename T>
      : Vector<T>,
	ConstantSizeAlgebraicCollection<T>
    {};
#else
    /// Concept ConstantSizeVector
    /**
        \par Refinement of:
	- Vector < T >
	- ConstantSizeAlgebraicCollection < T >
     */ 
    template <typename T>
    struct ConstantSizeVector
      : public Vector<T>,
	public ConstantSizeAlgebraicCollection<T>
    {};
#endif

    

/*@}*/ // end of group Concepts


} // namespace mtl

#endif // MTL_VECTOR_CONCEPTS_INCLUDE
