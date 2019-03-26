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

#ifndef MTL_DETAIL_RANGE_GENERATOR_INCLUDE
#define MTL_DETAIL_RANGE_GENERATOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/utility/complexity.hpp>
#include <boost/numeric/mtl/detail/base_cursor.hpp>
#include <boost/mpl/less.hpp>

namespace mtl { namespace traits { namespace detail {

    /// Range generator that traverses all elements of some densely stored collection 
    /** - Or contiguous parts of such collection
        - Works for matrices and vectors when derived from contiguous_memory_block
    **/
    template <typename Collection, typename Cursor, typename Complexity>
    struct dense_element_range_generator
    {
	typedef Complexity          complexity;
	typedef Cursor              type;
	static int const            level = 1;

	type begin(Collection const& collection)
	{
	    return collection.elements();
	}
	type end(Collection const& collection)
	{
	    return collection.elements() + collection.used_memory();
	}
    };

    /// Range generator that traverses all elements of some collection stored in strides
    template <typename Collection, typename Ref, typename Traversor>
    struct strided_element_range_generator
    {
	typedef complexity_classes::linear  complexity;
	typedef Traversor                   type;
	static int const                    level = 1;

	type begin(Ref& c)
	{
	    return type(c.address_data(), c.stride());
	}
	type end(Ref& c)
	{
		using mtl::size;
	    return type(c.address_data() + size(c) + c.stride(), c.stride());
	}
    };


    // Like above over all elements but in terms of offsets
    // Also with reference to collection in cursor
    template <typename Matrix, typename Cursor, typename Complexity>
    struct all_offsets_range_generator
    {
	typedef Complexity          complexity;
	typedef Cursor              type;
	static int const            level = 1;

	type begin(Matrix const& matrix) const
	{
	    return type(matrix, 0);
	}
	type end(Matrix const& matrix) const
	{
	    return type(matrix, matrix.nnz());
	}
    };
    

    // Cursor to some submatrix (e.g. row, column, block matrix, block row)
    // This cursor is intended to be used by range generators to iterate 
    // over subsets of the submatrix this cursor refers to.
    // For instance if this cursor refers to a row then a range 
    // can iterate over the elements in this row.
    // If this cursor refers to a block then a range can iterate over the rows in this block.
    // The level of a generated cursor must be of course at least one level less
    // The tag serves to dispatching between row and column cursors
    template <typename Matrix, typename Tag, int Level>
    struct sub_matrix_cursor
      : mtl::detail::base_cursor<std::size_t>
    {
	typedef sub_matrix_cursor                        self;
	typedef mtl::detail::base_cursor<std::size_t>    base;
	typedef Matrix                                   ref_type;
	static int const                                 level = Level;

	sub_matrix_cursor(std::size_t i, Matrix const& c)
	    : base(i), ref(c) 
	{}	

	self operator+(std::size_t offset) const
	{
	    return self(key + offset, ref);
	}

	// otherwise base_cursor returns an int and ranged for doesn't work
	// for getting the key of base_cursor use this->value()
	self operator*() const { return *this; }
	
	Matrix const& ref;
    };

    // Key for canonically referring to its elements with row and column index
    template <typename Matrix>
    struct matrix_element_key
    {
	typedef typename Collection<Matrix>::size_type           size_type;
	typedef matrix_element_key                               self;

	matrix_element_key(Matrix const& ref, std::size_t r, std::size_t c) : ref(ref)
	{
	    indices[0]= size_type(r); indices[1]= size_type(c);
	}

	bool operator==(const self& cc) const { return &ref == &cc.ref && indices[0] == cc.indices[0] && indices[1] == cc.indices[1]; }
	bool operator!=(const self& cc) const { return !(*this == cc); }

	size_type indices[2];
	Matrix const& ref;
    };

    // Cursor for canonically referring to its elements with row and column index
    // Increments row for pos==0 and column for pos==1
    // Referring operator returns matrix_element_key
    template <typename Matrix, int pos>
    struct matrix_element_cursor
    {
	typedef typename Collection<Matrix>::size_type           size_type;
	typedef matrix_element_cursor                            self;
	typedef matrix_element_key<Matrix>                       key_type;
	static int const            level = 2;
	
	matrix_element_cursor(Matrix const& ref, size_type r, size_type c) : ref(ref)
	{
	    indices[0]= r; indices[1]= c;
	}
	
	key_type operator*() const { return key_type(ref, indices[0], indices[1]); }

	self& operator++() { ++indices[pos]; return *this; }
	self operator++(int) { self tmp(*this); ++indices[pos]; return tmp; }
	template <typename T> self& operator+=(T n) { indices[pos] += size_type(n); return *this; }
	template <typename T> self& operator+(T n) const { self tmp = *this; tmp += n; return tmp; }

	self& operator--() { indices[pos]--; return *this; }
	self operator--(int) { self tmp(*this); indices[pos]--; return tmp; }
	template <typename T> self& operator-=(T n) { indices[pos] -= size_type(n); return *this; }
	template <typename T> self& operator-(T n) const { self tmp = *this; tmp -= n; return tmp; }

	bool operator==(const self& cc) const { return &ref == &cc.ref && indices[0] == cc.indices[0] && indices[1] == cc.indices[1]; }
	bool operator!=(const self& cc) const { return !(*this == cc); }

	size_type indices[2];
	Matrix const& ref;
    };


    template <typename Matrix, typename Complexity, int Level>
    struct all_rows_range_generator
    {
	typedef Complexity                                       complexity;
	static int const                                         level = Level;
	typedef Matrix                                           ref_type;
	typedef sub_matrix_cursor<Matrix, glas::tag::row, Level> type;
	typedef typename Collection<Matrix>::size_type           size_type;

	type begin(Matrix const& c) const
	{
	    return type(0, c); // return type(c.begin_row(), c); get rid of obsolete stuff
	}
	type end(Matrix const& c) const
	{
		using mtl::num_rows; using mtl::mat::num_rows;
	    return type(int(num_rows(c)), c); //return type(c.end_row(), c);
	}
	type lower_bound(Matrix const& c, size_type position) const
	{
		using mtl::num_rows;
	    return type(std::min(num_rows(c), position), c);
	}
    };

    template <typename Cursor>
    struct all_cols_in_row_range_generator
    {
	typedef complexity_classes::linear                       complexity;
	static int const                                         level = 1;
	typedef typename Cursor::ref_type                        ref_type;
	typedef typename Collection<ref_type>::size_type         size_type;


	typedef matrix_element_cursor<ref_type, 1>               type;

	type begin(Cursor const& c) const { return type(c.ref, c.value(), size_type(0)); }
	type end(Cursor const& c) const { using mtl::num_cols; return type(c.ref, c.value(), num_cols(c.ref)); }
	type lower_bound(Cursor const& c, size_type position) const
	{
		using mtl::num_cols;
	    return type(c.ref, c.value, std::min(num_cols(c.ref), position));
	}
    };


    template <typename Matrix, typename Complexity, int Level>
    struct all_cols_range_generator
    {
	typedef Complexity          complexity;
	static int const            level = Level;
	typedef sub_matrix_cursor<Matrix, glas::tag::col, Level> type;
	typedef typename Collection<Matrix>::size_type           size_type;

	type begin(Matrix const& c) const
	{
	    return type(0, c); // return type(c.begin_col(), c);
	}
	type end(Matrix const& c) const
	{
		using mtl::num_cols;
	    return type(num_cols(c), c); // return type(c.end_col(), c);
	}
	type lower_bound(Matrix const& c, size_type position) const
	{
		using mtl::num_cols;
	    return type(std::min(num_cols(c), position), c);
	}
    };

    template <typename Cursor>
    struct all_rows_in_col_range_generator
    {
	typedef complexity_classes::linear                       complexity;
	static int const            level = 1;
	typedef typename Cursor::ref_type                        ref_type;
	typedef typename Collection<ref_type>::size_type         size_type;


	typedef matrix_element_cursor<ref_type, 0> type;

	type begin(Cursor const& c) const { return type(c.ref, 0, c.value()); }
	type end(Cursor const& c) const { using mtl::num_rows; return type(c.ref, num_rows(c.ref), c.value()); }
	type lower_bound(Cursor const& c, size_type position) const
	{
		using mtl::num_rows;
	    return type(c.ref, std::min(num_rows(c.ref), position), c.value());
	}
    };


    // Use RangeGenerator for Collection by applying to .ref
    template <typename Coll, typename RangeGenerator>
    struct referred_range_generator
    {
	typedef typename RangeGenerator::complexity complexity;
	static int const                            level = RangeGenerator::level;
	typedef typename RangeGenerator::type       type;
	typedef typename Collection<Coll>::size_type  size_type;
	
	type begin(const Coll& c)
	{
	    return RangeGenerator().begin(c.ref);
	}

	type end(const Coll& c)
	{
	    return RangeGenerator().end(c.ref);
	}
	type lower_bound(const Coll& c, size_type position)
	{
	    return RangeGenerator().lower_bound(c.ref, position);
	}
    };

} // namespace detail

    namespace range {

	template <typename Range1, typename Range2>
	struct min
	    : public boost::mpl::if_< 
	            boost::mpl::less<
	                   typename Range1::complexity, 
	                   typename Range2::complexity>
	          , Range1
	          , Range2
	      >
	{};
    }

}} // namespace mtl::traits

#endif // MTL_DETAIL_RANGE_GENERATOR_INCLUDE
