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

#ifndef MTL_MATRIX_CONCEPT_INCLUDE
#define MTL_MATRIX_CONCEPT_INCLUDE

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
    concept Matrix<typename T>
      : AlgebraicCollection<T>
    {
	const_reference T::operator() (size_type row, size_type col) const;

	size_type nnz(T);

	// A[r][c] equivalent to A(r, c)
    };
#else
    /// Concept Matrix
    /**
        \par Refinement of:
	- AlgebraicCollection < T >
	\par Notation:
	- X is a type that models Matrix
	- A is an object of type X
	- r, c are objects of size_type
	\par Valid expressions:
	- Element access: \n A(r, c) \n Return Type: const_reference 
	  \n Semantics: Element in row \p r and column \p c
	- Element access: \n A[r][c] \n Equivalent to A(r, c)
	\par Models:
	- dense2D
	- morton_dense
	- compressed2D
	\note
	-# The access via A[r][c] is supposed to be implemented by means of A(r, c) (typically via CRTP and proxies).
	  If it would become (extremely) important to support 2D C arrays, it might be necessary to drop the requirement
	  of element access by A(r, c).
	-# The name const_reference does not imply that the return type is necessarily referrable. For instance compressed2D
	   returns value_type.
     */ 
    template <typename T>
    struct Matrix
	: public AlgebraicCollection<T>
    {
	/// Element access
	const_reference T::operator() (size_type row, size_type col) const;
    };
#endif

    

#ifdef __GXX_CONCEPTS__
    concept MatrixInserter<typename T>
    {
	typename matrix_type;
	// typename T::matrix_type;

	requires Matrix<matrix_type>;
	
	typename proxy_type;
	proxy_type operator() (Matrix<matrix_type>::size_type row, Matrix<matrix_type>::size_type col);
	
	T operator<< (proxy_type, Matrix<matrix_type>::value_type>);
    };
#else
    /// Concept MatrixInserter: classes that enable efficient insertion into matrices, esp. compressed sparse.
    /** 
	Used to fill non-mutable matrices like compressed2D. Matrix inserters might be parametrizable with
	update functor. This allow to perform different operations when entry already exist, e.g. overwriting,
	incrementing, minimum, ... The most important updates are certainly overwrite and increment (add).

	\par Associated types
	- matrix_type

	\par Requires:
	- Matrix<matrix_type>
	
	\par Notation:
	- X is a type that models MatrixInserter
	- A is an object of type X
	- r, c are objects of type Matrix<matrix_type>::size_type
	- v is an object of type Matrix<matrix_type>::value_type

	\par Valid expressions:
	- Insertion with shift operator: \n
	   A(r, c) << v \n
	   Return type: T
	\par Models:
	- mtl::mat::inserter < T >
	\note
	-# Used in concept InsertableMatrix
     */
    template <typename T>
    struct MatrixInserter
    {
	/// Type  of matrix into which is inserted
	typedef associated_type matrix_type;

	/// Return type of element access; only proxy
	typedef associated_type  proxy_type;
	/// Element access; returns a proxy that handles insertion
	proxy_type operator() (Matrix<matrix_type>::size_type row, Matrix<matrix_type>::size_type col);
    };
#endif

#ifdef __GXX_CONCEPTS__
    concept InsertableMatrix<typename T>
      : Matrix<T>
    {
	requires MatrixInserter<mtl::mat::inserter<T> >;
    };
#else
    /// Concept InsertableMatrix: %matrix that can be filled by means of inserter
    /** 
	\par Requires:
	- MatrixInserter < mtl::mat::inserter< T > >
	\par Models:
	- dense2D
	- morton_dense
	- compressed2D
	\note
	-# All matrices in MTL model this concept in order and all future matrices are supposed to.
    */
    template <typename T>
    struct InsertableMatrix
      : Matrix < T >
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept MutableMatrix<typename T>
      : Matrix<T>,
	MutableCollection<T>
    {
	reference T::operator() (size_type row, size_type col);

	// A[r][c] equivalent to A(r, c)
    };
#else
    /// Concept MutableMatrix
    /**
        \par Refinement of:
	- Matrix < T >
	- MutableCollection < T >
	\par Notation:
	- X is a type that models MutableMatrix
	- A is an object of type X
	- r, c are objects of size_type
	\par Valid expressions:
	- Element access: \n A(r, c) \n Return Type: reference 
	  \n Semantics: Element in row \p r and column \p c
	- Element access: \n A[r][c] \n Equivalent to A(r, c)
	\par Models:
	- dense2D
	- morton_dense
	\note
	-# The access via A[r][c] is supposed to be implemented by means of A(r, c) (typically via CRTP and proxies).
	  If it would become (extremely) important to support 2D C arrays, it might be necessary to drop the requirement
	  of element access by A(r, c).
     */ 
    template <typename T>
    struct MutableMatrix
	: public Matrix<T>,
	  public MutableCollection<T>
    {
	/// Element access (in addition to const access)
	reference T::operator() (size_type row, size_type col);
    };
#endif

    
#ifdef __GXX_CONCEPTS__
    concept ConstantSizeMatrix<typename T>
      : Matrix<T>,
	ConstantSizeAlgebraicCollection<T>
    {};
#else
    /// Concept ConstantSizeMatrix
    /**
        \par Refinement of:
	- Matrix < T >
	- ConstantSizeAlgebraicCollection < T >
     */ 
    template <typename T>
    struct ConstantSizeMatrix
      : public Matrix<T>,
	public ConstantSizeAlgebraicCollection<T>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept ResizeableMatrix<typename T>
      : Matrix<T>
    {
	void T::resize(size_type r, size_type c);
    };
#else
    /// Concept ResizeableMatrix
    /**
        \par Refinement of:
	- Matrix < T >
     */ 
    template <typename T>
    struct ResizeableMatrix
      : public Matrix<T>
    {
	/// Resize function
	/** If new memory should be allocated only if total size changes */
	void resize(size_type r, size_type c);
    };
#endif


#ifdef __GXX_CONCEPTS__
    concept RowTraversableMatrix<typename M>
      : Matrix<M>,
        TraversableCollection<mtl::tag::row, M> 
    {};
#else
    /// Concept RowTraversableMatrix: provides begin and end cursor to traverse rows
    /**
        \par Refinement of:
	- Matrix < M >
	- TraversableCollection <mtl::tag::row, M> 
     */ 
    template <typename M>
    struct RowTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<mtl::tag::row, M>
    {};
#endif

    
#ifdef __GXX_CONCEPTS__
    concept ColumnTraversableMatrix<typename M>
      : Matrix<M>,
        TraversableCollection<mtl::tag::col, M> 
    {};
#else
    /// Concept ColumnTraversableMatrix: provides begin and end cursor to traverse columns
    /**
        \par Refinement of:
	- Matrix < M >
	- TraversableCollection <mtl::tag::col, M> 
     */ 
    template <typename M>
    struct ColumnTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<mtl::tag::col, M>
    {};
#endif


#ifdef __GXX_CONCEPTS__
    concept MajorTraversableMatrix<typename M>
      : Matrix<M>,
        TraversableCollection<mtl::tag::major, M> 
    {};
#else
    /// Concept MajorTraversableMatrix: traversable on major dimension
    /**
        Concept for matrices that are traversable along the major dimension, i.e.
	traversing the rows of a row-major matrix and the columns of a column-major matrices.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
	- TraversableCollection <mtl::tag::major, M> 
	\note
	-# This traversal corresponds to the iterator design in MTL 2.
     */ 
    template <typename M>
    struct MajorTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<mtl::tag::major, M>
    {};
#endif
    

#ifdef __GXX_CONCEPTS__
    concept MinorTraversableMatrix<typename M>
      : Matrix<M>,
        TraversableCollection<mtl::tag::minor, M> 
    {};
#else
    /// Concept MinorTraversableMatrix: traversable on minor dimension
    /**
        Concept for matrices that are traversable along the minor dimension, i.e.
	traversing the columns of a row-major matrix and the rows of a column-major matrices.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
	- TraversableCollection <mtl::tag::minor, M> 
	\note
	-# This traversal corresponds to the iterator design in MTL 2.
     */ 
    template <typename M>
    struct MinorTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<mtl::tag::minor, M>
    {};
#endif
    

#ifdef __GXX_CONCEPTS__
    concept AllTraversableMatrix<typename M>
      : Matrix<M>,
        TraversableCollection<mtl::tag::all, M> 
    {};
#else
    /// Concept AllTraversableMatrix: provides traversion over all elements
    /**
        All elements of a matrix are traversed, including structural zeros. Can be used, e.g.,
	for printing.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
	- TraversableCollection <mtl::tag::all, M> 
	\note
	-# For dense matrices the concept is equivalent to NonZeroTraversableMatrix.
     */ 
    template <typename M>
    struct AllTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<mtl::tag::all, M>
    {};
#endif
    

#ifdef __GXX_CONCEPTS__
    concept NonZeroTraversableMatrix<typename M>
      : Matrix<M>,
        TraversableCollection<mtl::tag::nz, M> 
    {};
#else
    /// Concept NonZeroTraversableMatrix: provides traversion over all structural non-zeros
    /**
        All structural non-zero elements of a matrix are traversed. Can be used, e.g.,
	for copying.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
	- TraversableCollection <mtl::tag::all, M> 
	\note
	-# For dense matrices the concept is equivalent to AllTraversableMatrix.
     */ 
    template <typename M>
    struct AllTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<mtl::tag::all, M>
    {};
#endif
    

#ifdef __GXX_CONCEPTS__
    concept AllTraversableSubMatrix<typename Tag, typename M>
      : Matrix<M>,
        TraversableCollection<Tag, M>,
	TraversableCollection<mtl::tag::all, TraversableCollection<Tag, M>::result_type>
    {};
#else
    /// Concept AllTraversableSubMatrix: provides traversion of rows, columns of matrices
    /**
        All elements of a row or a column, according to the Tag, are traversed.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
        - TraversableCollection<Tag, M>,
	- TraversableCollection<mtl::tag::all, TraversableCollection<Tag, M>::result_type>
     */ 
    template <typename Tag, typename M>
    struct AllTraversableMatrix
      : public Matrix<M>,
        public TraversableCollection<Tag, M>,
	public TraversableCollection<mtl::tag::all, TraversableCollection<Tag, M>::result_type>
    {};
#endif
    
    

#ifdef __GXX_CONCEPTS__
    concept NonZeroTraversableSubMatrix<typename Tag, typename M>
      : Matrix<M>,
        TraversableCollection<Tag, M>,
	TraversableCollection<mtl::tag::nz, TraversableCollection<Tag, M>::result_type>
    {};
#else
    /// Concept NonZeroTraversableSubMatrix: provides traversion of non-zero in rows or columns of matrices
    /**
        All structural non-zero elements of a row or a column, according to the Tag, are traversed.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
        - TraversableCollection<Tag, M>,
	- TraversableCollection<mtl::tag::nz, TraversableCollection<Tag, M>::result_type>
     */ 
    template <typename Tag, typename M>
    struct NonZeroTraversableSubMatrix
      : public Matrix<M>,
        public TraversableCollection<Tag, M>,
	public TraversableCollection<mtl::tag::nz, TraversableCollection<Tag, M>::result_type>
    {};
#endif
    

#ifdef __GXX_CONCEPTS__
    concept IteratableSubMatrix<typename Tag, typename ITag, typename M>
      : Matrix<M>,
        TraversableCollection<Tag, M>,
	TraversableCollection<ITag, TraversableCollection<Tag, M>::result_type>
    {};
#else
    /// Concept IteratableSubMatrix: provides iteration over elements within rows or columns of matrices
    /**
        This concepts actually combines four sub-concepts. The iteration can be either performed over
	all elements or only over structural non-zero elements whereby the iterator can be a const-iterator
	or a mutable iterator. These four combinations are specified by the tags mtl::tag::iter::all, 
	mtl::tag::iter::nz, mtl::tag::const_iter::all,  and
	mtl::tag::const_iter::nz for ITag. The template parameter Tag can be mtl::tag::major or mtl::tag::column.
	The cursors begin and end are provided.
        \par Refinement of:
	- Matrix < M >
        - TraversableCollection<Tag, M>,
	- TraversableCollection<ITag, TraversableCollection<Tag, M>::result_type>
     */ 
    template <typename Tag, typename ITag, typename M>
    struct IteratableSubMatrix
      : public Matrix<M>,
        public TraversableCollection<Tag, M>,
	public TraversableCollection<ITag, TraversableCollection<Tag, M>::result_type>
    {};
#endif






/*@}*/ // end of group Concepts

} // namespace mtl

#endif // MTL_MATRIX_CONCEPT_INCLUDE
