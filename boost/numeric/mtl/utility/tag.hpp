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

#ifndef MTL_TAG_INCLUDE
#define MTL_TAG_INCLUDE

#include <cstddef>
#include <boost/numeric/mtl/utility/glas_tag.hpp>


namespace mtl { namespace tag {

/** @defgroup Tags Tags for concept-free dispatching
 *  @{
 */

// For internal use (e.g., to invalidate range generators)
// Is this the same as bottom?
struct unsupported {};

// Name says it
struct dummy3 {};
struct dummy4 {};

/// Tag for all types
struct universe {};

/// Tag used to flatten categories
/** The virtual derivation causes perceivable run-time overhead that can be avoided with this struct using traits::flatcat1 and such. **/
template <typename T> struct flat : universe {};

// Tag for any scalar value
/** At the moment default for unknown types (will be precised later) */
struct scalar : virtual universe {};

/// For non-MTL types with category not explicitly defined
/** At the moment all treated as scalars (will be precised later) */
struct unknown : virtual scalar {};

/// Tag for intermediate objects that require explicit evaluation 
struct unevaluated : virtual universe {};

/// Any collection, i.e. vectors, matrices or higher-dimensional tensor
struct collection : virtual universe {};

/// Tag for any MTL vector (and user-declared MTL vectors)
struct vector : virtual collection {};

/// Tag for references to vector
/** For instance to not access memory directly but use functions, e.g. in set_to_zero. **/ 
struct vector_ref : virtual vector {};

/// Tag for any MTL column vector (and user-declared MTL vectors)
struct col_vector : virtual vector {};

/// Tag for any MTL row vector (and user-declared MTL vectors)
struct row_vector : virtual vector {};

/// Tag for any MTL matrix (and user-declared MTL matrices)
struct matrix : virtual collection {};

/// Tag for any dense collection
struct dense : virtual universe {};
    
/// Tag for vectors with one-dimensional memory addressing
/** offet v_i is x*i for some x */
struct has_1D_layout : virtual dense {};
    
/// Tag for matrices with two-dimensional memory addressing
/** offet a_ij is x*i + y*j for some x and y */
struct has_2D_layout : virtual dense {};
    
/// Tag for any sparse collection
struct sparse : virtual universe {};
    
// for distinction between dense and sparse matrices
struct dense_matrix : virtual dense, virtual matrix {};
struct sparse_matrix : virtual sparse, virtual matrix {};

/// Tag for collections where values are stored contigously in memory
struct contiguous_memory : virtual universe {};

/// Tag for dense and contiguous collections
/** Only short cut */
struct contiguous_dense : virtual dense, virtual contiguous_memory {};

/// Collection with iterator
struct has_iterator : virtual universe {};

/// Collection with random-access iterator
struct has_ra_iterator : virtual has_iterator {};

/// Collection with fast random-access iterator
/** Meaning: unrolling is probably beneficial. Counter-example: Morton-ordered matrices have
    random access but this is so slow that regular traversal is favorable */
struct has_fast_ra_iterator : virtual has_ra_iterator {};

/// Collection with cursor
struct has_cursor : virtual universe {};

/// Collection with random-access cursor
struct has_ra_cursor : virtual has_cursor {};

/// Collection with fast random-access cursor
/** Meaning: unrolling is probably beneficial. Counter-example: Morton-ordered matrices have
    random access but this is so slow that regular traversal is favorable */
struct has_fast_ra_cursor : virtual has_ra_cursor {};

/// Tag for matrices with sub_matrix function exist and doesn't say for which ranges it is defined
struct has_sub_matrix : virtual universe {};

/// Sub-divisible into quadrants, i.e. arbitrary sub-matrices not necessarily supported but recursion works
//  more explanation needed
struct qsub_divisible : virtual has_sub_matrix {};

/// Tag for sub-divisible matrix, i.e. sub_matrix works 
struct sub_divisible : virtual qsub_divisible {};

/// Tag for dense row vector in the category lattice
struct dense_row_vector
  : virtual row_vector, virtual contiguous_dense, 
    virtual has_fast_ra_iterator, virtual has_fast_ra_cursor, virtual has_1D_layout
{};

/// Tag for dense column vector in the category lattice
struct dense_col_vector
  : virtual col_vector, virtual contiguous_dense, 
    virtual has_fast_ra_iterator, virtual has_fast_ra_cursor, virtual has_1D_layout
{};

/// Tag for strided row vector in the category lattice
struct strided_row_vector
  : virtual row_vector, virtual vector_ref,
    virtual has_fast_ra_iterator, virtual has_fast_ra_cursor, virtual has_1D_layout
{};

/// Tag for strided column vector in the category lattice
struct strided_col_vector
  : virtual col_vector, virtual vector_ref, 
    virtual has_fast_ra_iterator, virtual has_fast_ra_cursor, virtual has_1D_layout
{};

/// Tag for sparse row vector in the category lattice
struct sparse_row_vector
  : virtual row_vector, virtual sparse
{};

/// Tag for sparse column vector in the category lattice
struct sparse_col_vector
  : virtual col_vector, virtual sparse
{};

/// Tag to handle std::vector in the category lattice
struct std_vector
  : virtual vector, virtual contiguous_dense, virtual has_1D_layout
{};

/// Tag for a view on a (regular) dense matrix in the category lattice
/** The map perform address computation and has therefore no 2D-layout.
    It is also not (yet) assumed that the view provides iterators. */
struct dense2D_view 
  : virtual matrix, virtual contiguous_dense, virtual has_fast_ra_cursor 
    // , virtual sub_divisible // is currently not sub-divisible
{};

/// Tag for (regular) dense matrix in the category lattice
struct dense2D 
  : virtual dense2D_view, virtual has_fast_ra_iterator, virtual has_2D_layout, virtual sub_divisible
{};

struct implicit_dense
  : virtual matrix, virtual dense, virtual has_fast_ra_cursor
{};

/// Tag for a view on a Morton-order matrix in the category lattice
/** It is not (yet) assumed that the view provides iterators. */
struct morton_view 
  : virtual matrix, virtual contiguous_dense,  
    virtual has_ra_cursor // , virtual qsub_divisible // is currently not sub-divisible
{};


/// Tag for Morton-order matrix in the category lattice
struct morton_dense 
  : virtual morton_view, virtual has_ra_iterator, virtual qsub_divisible
{};

/// Tag for a view on a compressed matrix in the category lattice
/** It is not (yet) assumed that the view provides iterators. */
struct compressed2D_view
  : virtual matrix, virtual sparse, virtual has_cursor
{};

/// Tag for compressed matrix in the category lattice
struct compressed2D 
  : virtual compressed2D_view, virtual has_iterator
{};

/// Tag for multi_vector
// Maybe splitting later into sparse and dense form
struct multi_vector
  : virtual matrix, virtual dense
{};

/// Tag for transposed multi_vector
// Maybe splitting later into sparse and dense form
struct transposed_multi_vector
  : virtual matrix, virtual dense
{};

/// Tag for transposed multi_vector
// Maybe splitting later into sparse and dense form
struct hermitian_multi_vector
  : virtual matrix, virtual dense
{};

/// Tag for ell_matrix (preliminary)
struct ell_matrix
  : sparse_matrix
{};

/// Tag for element structure matrix
struct element_structure
  : sparse_matrix
{};

/// Tag for element structure matrix
struct sparse_banded_matrix
  : sparse_matrix
{};

/// Tag for mat_cvec_multiplier
struct mat_cvec_multiplier
  : col_vector
{};

/// Tag for implicit dense matrices


/// Tag for bottom of the category lattice
/** Only for completeness; probably not needed in practice. */
struct bottom
  : virtual compressed2D, virtual morton_dense, virtual dense2D, 
    virtual dense_col_vector, virtual dense_row_vector, virtual unknown,
    virtual multi_vector
{};

template <typename Tag1, typename Tag2, typename Tag3= dummy3, typename Tag4= dummy4>
struct join
  : virtual Tag1, virtual Tag2, virtual Tag3, virtual Tag4
{};


// =====================
// Types for orientation
// =====================


/// Characterizes row-major orientation in matrices and row vector in 1D
struct row_major {};

/// Characterizes column-major orientation in matrices and column vector in 1D
struct col_major {};

/// Synonym for col_major
typedef col_major column_major;

/// Common base for diagonal tags
struct universe_diagonal {};

/// Tag indicating that diagonal is stored regularly
struct regular_diagonal : universe_diagonal {};

/// Tag indicating that diagonal contains unit elements
struct unit_diagonal : universe_diagonal {};

/// Tag indicating that diagonal entries are stored as inverses
/** Storing value in different ways can be faster in several algorithms.
    By the time of this writing it is experimental and only used
    in upper_trisolve and lower_trisolve. **/
struct inverse_diagonal : universe_diagonal {};




/*@}*/ // end of group Tags

} // namespace mtl::tag

/** @addtogroup Tags
 *  @{
 */

// Import in mtl namespace
using tag::row_major;
using tag::col_major;
using tag::column_major;
using tag::sparse;
using tag::dense;

#ifdef MTL_DOX_ONLY
// using is not documented therefor the redeclaration here

/// Characterizes row-major orientation in matrices and row vector in 1D, import of tag::row_major
struct row_major {};

/// Characterizes column-major orientation in matrices and column vector in 1D, import of tag::col_major
struct col_major {};

/// Synonym for col_major, import of tag::column_major
struct column_major {};

/// Tag for any sparse collection, import of tag::sparse
struct sparse {};

/// Tag for any dense collection, import of tag::dense
struct dense {};

#endif


/// Tag for any banded collection
struct banded {};

/// Tag for any compressed collection
struct compressed {};

/// Tag for coordinate matrix
struct coordinate {};

/// Tag for Ellpack matrices
struct ellpack {};

/// Tag synonym for Ellpack matrices
typedef ellpack ell;

/// Tag for any morton-orded collection
struct morton {};

/// Tag for any symmetric collection
struct symmetric {};

/// Tag for any anti-symmetric collection
struct anti_symmetric {};

/// Tag for any self-adjoint, i.e. Hermitian collection
struct self_adjoint {};

/// Tag for specifying size type in type generators
template <typename T>
struct as_size_type 
{
    typedef T type;
};

/// Tag for any collection stored on stack
struct on_stack {};

/// Tag for any collection stored on heap
struct on_heap {};

#ifdef MTL_WITH_VARIADIC_TEMPLATE

/// To define compile-time dimension in type generators
template <std::size_t ...Values>
struct dim {};

/// To define mask for Morton-order 
template <std::size_t ...Values>
struct mask {};

#endif

/*@}*/ // end of group Tags

// =====================
// Tags for traversals
// Import some from GLAS
// =====================

namespace tag {

/** @addtogroup Tags
 *  @{
 */

    /// Tag for cursor traversal of non-zero elements of a collection
    /** Also used for elements within rows and columns */
    using glas::tag::nz;

    /// Tag for cursor traversal of all elements of a collection
    /** Also used for elements within rows and columns */
    using glas::tag::all;

    /// Tag for traversal of all rows in a matrix
    using glas::tag::row;
    /// Tag for traversal of all columns in a matrix
    using glas::tag::col;

    /// Tag for traversal of a matrix's major dimension
    /** Is equivalent to glas::tag::row for row-major matrices and 
	glas::tag::col for column-major matrices */
    using glas::tag::major;

    /// Tag for traversal of a matrix's minor dimension
    /** Is equivalent to glas::tag::row for row-major matrices and 
	glas::tag::col for column-major matrices */
    using glas::tag::minor;

    // To define iterators over matrices or rows/cols of it, vectors

    namespace iter {

	/// Tag for iterator traversal of non-zero elements of a collection
	/** Also used for elements within rows and columns */
	struct nz {};

	/// Tag for iterator traversal of all elements of a collection
	/** Also used for elements within rows and columns */
	struct all {};

    } // namespace mtl::tag::iter

    // Same with const iterators

    namespace const_iter {

	/// Tag for const-iterator traversal of non-zero elements of a collection
	/** Also used for elements within rows and columns */
	struct nz {};

	/// Tag for const-iterator traversal of all elements of a collection
	/** Also used for elements within rows and columns */
	struct all {};

    } // namespace mtl::tag::const_iter

/*@}*/ // end of group Tags

} // namespace mtl::tag


} // namespace mtl

#endif // MTL_TAG_INCLUDE
