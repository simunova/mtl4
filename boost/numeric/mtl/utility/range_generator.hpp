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

#ifndef MTL_RANGE_GENERATOR_INCLUDE
#define MTL_RANGE_GENERATOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/detail/range_generator.hpp>
#include <boost/numeric/mtl/utility/complexity.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

namespace mtl {

namespace traits
{
    /// Functor for generating begin and end cursors over a collection
    /** Thus functor must contain begin and end member functions which are used by free functions
        The cursor type must be defined as 'typedef xxx type;'
        complexity characterizes the run time of a traveral, cf. utility/complexity.hpp
          complexity can be used to dispatch between different traversals depending on algorithm
          and on collection (corr. subset represented by cursor)
        level indicates the maximal level of nesting
          - level 0: there is no traversal of this Tag for this collection
          - level 1: cursor refers directly to elements
          - level 2: cursor iterates over sets of elements or is only an intermediate cursor
                     cursor e.g. over rows, 
                     its generated ranges are level 1 and iterate over elements
          - level 3: cursor over sets of sets of elements,
                     its generated ranges are level 2 or 1 depending on the tag used on the cursor
          - level 4: for instance blocked matrix -> level 4: block rows -> level 3: block elements
                       -> level 2: regular rows -> level 1: matrix elements
                     if an element cursor range was generated from the block element then the nesting 
                       would be only 3 (since the last two levels collapse)
        Cursors of level > 1 represent subsets of a collection and thus, it is only logical that
        there must be range generators for these subset, which are applied on the cursor. 
    **/
    template <typename Tag, typename Collection>
    struct range_generator
    {
	typedef complexity_classes::infinite  complexity;
	static int const             level = 0;
	typedef Tag                  tag;
	
	// MTL_STATIC_ASSERT(false, "range generator not implemented for these arguments (e.g. wrong combination).");

	// specializations must contain the following members
	// typedef xxx               type;
	// type begin() { ... }
	// type end()   { ... }
    };
} // namespace traits



/// Returns begin cursor over the Collection or a subset of the Collection
/** Form of traversal depends on Tag, cf utility/glas_tag.hpp
    On nested traversals, cursors of level > 1 must provide at least one range generator 
    \sa end() **/
template <class Tag, class Collection>
typename traits::range_generator<Tag, Collection>::type 
inline begin(Collection const& c)
{
    typedef traits::range_generator<Tag, Collection> gen_type;
    MTL_STATIC_ASSERT(gen_type::level != 0, 
		      "Template arguments not supported, probably traversal with unsupported tag combination.");
    return gen_type().begin(c);
}

template <class Tag, class Collection>
typename traits::range_generator<Tag, Collection>::type 
inline begin(Collection& c)
{
    typedef traits::range_generator<Tag, Collection> gen_type;
    MTL_STATIC_ASSERT(gen_type::level != 0, 
		      "Template arguments not supported, probably traversal with unsupported tag combination.");
    return gen_type().begin(c);
}

/// Corresponding end cursor
/** \sa begin() **/
template <class Tag, class Collection>
typename traits::range_generator<Tag, Collection>::type 
inline end(Collection const& c)
{
    return traits::range_generator<Tag, Collection>().end(c);
}

template <class Tag, class Collection>
typename traits::range_generator<Tag, Collection>::type 
inline end(Collection& c)
{
    return traits::range_generator<Tag, Collection>().end(c);
}

/// Cursor to an element with index equal or larger than \p position in a one-dimensional traversion.
/** This function is only defined where Tag represents an obvious one-dimensional traversion
    of Collection allowing for an interpretation of position. 
    Examples are tag::row for a dense matrix or row-major compressed matrix.
    tag::all regarding an entire matrix (i.e. going over all entries of a matrix)
    does not characterize a one-dimensional traversion
    so that \p position has no unique meaning.
    Traversing all entries of a matrix row is one-dimensional and \p position is understood as looking for
    a column index, i.e. lower_bound<tag::all>(row_cursor, 7) returns a cursor to a matrix element whose row
    is the one of row_cursor and whose column is 7.
    Likewise traversing all non-zeros of a row will return a cursor that points to an entry with according
    row index and column index at least 7.
    If the searched index could not be found in the one-dimensional collection the returned cursor will be 
    identical to the one returned by the end function.
**/
template <class Tag, class Coll>
typename traits::range_generator<Tag, Coll>::type 
inline lower_bound(Coll const& c, typename Collection<Coll>::size_type position)
{
    return traits::range_generator<Tag, Coll>().lower_bound(c, position);
}

template <class Tag, class Coll>
typename traits::range_generator<Tag, Coll>::type 
inline lower_bound(Coll& c, typename Collection<Coll>::size_type position)
{
    return traits::range_generator<Tag, Coll>().lower_bound(c, position);
}


namespace traits {

    // Dispatch between row and column-major traversal
    template <typename Matrix>
    struct range_generator<tag::major, Matrix>
	: public range_generator<
	    typename boost::mpl::if_<
	       boost::is_same<typename OrientedCollection<Matrix>::orientation, row_major>
	     , ::mtl::tag::row
             , ::mtl::tag::col
            >::type, 
          Matrix>
    {};

    // Dispatch between row and column-major traversal
    template <typename Matrix>
    struct range_generator<tag::minor, Matrix>
	: public range_generator<
	    typename boost::mpl::if_<
	       boost::is_same<typename OrientedCollection<Matrix>::orientation, row_major>
             , ::mtl::tag::col
	     , ::mtl::tag::row
            >::type, 
          Matrix>
    {};

    //=== Range generators for mat::indirect ====================

    // Range generator over all rows
    template <typename Matrix>
    struct range_generator<glas::tag::row, mtl::mat::indirect<Matrix> >
      : detail::all_rows_range_generator<mtl::mat::indirect<Matrix>, complexity_classes::linear, 2>
    {};

    // Range generator over all entries within a row
    template <typename Matrix>
    struct range_generator<glas::tag::all, 
			   mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::row, 2> >
    : detail::all_cols_in_row_range_generator<
	mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::row, 2> >
    {};

    // Range generator over all non-zero entries within a row, same as all entries
    template <typename Matrix>
    struct range_generator<glas::tag::nz, 
			   mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::row, 2> >
    : detail::all_cols_in_row_range_generator<
	mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::row, 2> >
    {};

    // - same for columns first

    // Range generator over all colums
    template <typename Matrix>
    struct range_generator<glas::tag::col, mtl::mat::indirect<Matrix> >
      : detail::all_cols_range_generator<mtl::mat::indirect<Matrix>, complexity_classes::linear, 2>
    {};

    // Range generator over all entries within a column
    template <typename Matrix>
    struct range_generator<glas::tag::all, 
			   mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::col, 2> >
    : detail::all_rows_in_col_range_generator<
	mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::col, 2> >
    {};

    // Range generator over all non-zero entries within a column, same as all entries
    template <typename Matrix>
    struct range_generator<glas::tag::nz, 
			   mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::col, 2> >
    : detail::all_rows_in_col_range_generator<
	mtl::traits::detail::sub_matrix_cursor<mtl::mat::indirect<Matrix>, glas::tag::col, 2> >
    {};

    // Take major same as row
    // Not necessarily best choice but shouldn't matter here
    template <typename Matrix>
    struct range_generator<glas::tag::major, mtl::mat::indirect<Matrix> >
      : detail::all_rows_range_generator<mtl::mat::indirect<Matrix>, complexity_classes::linear, 2>
    {};

}


} // namespace mtl

#endif // MTL_RANGE_GENERATOR_INCLUDE
