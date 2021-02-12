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

#ifndef MTL_DENSE2D_INCLUDE
#define MTL_DENSE2D_INCLUDE


#include <algorithm>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>
#include <boost/numeric/mtl/matrix/base_sub_matrix.hpp>
#include <boost/numeric/mtl/matrix/all_mat_expr.hpp>
#include <boost/numeric/mtl/matrix/operators.hpp>
#include <boost/numeric/mtl/detail/contiguous_memory_block.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/compute_factors.hpp>
#include <boost/numeric/mtl/operation/clone.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/utility/common_include.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/mtl/utility/dense_el_cursor.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/strided_dense_el_cursor.hpp>
#include <boost/numeric/mtl/utility/strided_dense_el_iterator.hpp>
#include <boost/numeric/mtl/utility/transposed_orientation.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

#ifdef MTL_WITH_INITLIST
# include <initializer_list>
#endif
#include <algorithm>

// Forward declaration (for friend declaration)
namespace mtl { namespace traits { namespace detail {
    template <typename, typename, bool> struct dense2D_iterator_range_generator;
}}}

namespace mtl { namespace mat {


using std::size_t;

// Forward declarations
template <typename Value, typename Parameters> class dense2D;
class dense2D_indexer;

// Helper type
struct dense2D_sub_ctor {};

// Indexing for dense matrices
class dense2D_indexer 
{
    // helpers for public functions
    size_t offset(size_t ldim, size_t r, size_t c, row_major) const 
    {
	return r * ldim + c; 
    }
    size_t offset(size_t ldim, size_t r, size_t c, col_major) const 
    {
	return c * ldim + r; 
    }
    
    size_t row(size_t offset, size_t ldim, row_major) const 
    {
	return offset / ldim; 
    }
    size_t row(size_t offset, size_t ldim, col_major) const 
    {
	return offset % ldim;
    }
    
    size_t col(size_t offset, size_t ldim, row_major) const 
    {
	return offset % ldim;
    }
    size_t col(size_t offset, size_t ldim, col_major) const 
    {
	return offset / ldim; 
    }

 public:
    template <typename Value, class Parameters>
    size_t operator() (const dense2D<Value, Parameters>& ma, size_t r, size_t c) const
    {
	typedef dense2D<Value, Parameters> matrix_type;
	// convert into c indices
	typename matrix_type::index_type my_index;
	size_t my_r= index::change_from(my_index, r);
	size_t my_c= index::change_from(my_index, c);
	return offset(ma.ldim, my_r, my_c, typename matrix_type::orientation());
    }

    template <typename Value, class Parameters>
    size_t row(const dense2D<Value, Parameters>& ma, 
	       typename dense2D<Value, Parameters>::key_type key) const
    {
	typedef dense2D<Value, Parameters> matrix_type;
	// row with c-index for my orientation
	size_t r= row(ma.offset(key), ma.ldim, typename matrix_type::orientation());
	return index::change_to(typename matrix_type::index_type(), r);
    }

    template <typename Value, class Parameters>
    size_t col(const dense2D<Value, Parameters>& ma, 
	       typename dense2D<Value, Parameters>::key_type key) const 
    {
	typedef dense2D<Value, Parameters> matrix_type;
	// column with c-index for my orientation
	size_t c= col(ma.offset(key), ma.ldim, typename matrix_type::orientation());
	return index::change_to(typename matrix_type::index_type(), c);
    }
    template <typename, typename> friend class dense2D;
}; // dense2D_indexer


namespace detail 
{
    
    // Compute required memory
    // Enabling mechanism to make sure that computation is valid
    template <typename Parameters, bool Enable>
    struct dense2D_array_size {
	static std::size_t const value= 0;
    };

    template <typename Parameters>
    struct dense2D_array_size<Parameters, true>
    {
	typedef typename Parameters::dimensions   dimensions;
	MTL_STATIC_ASSERT((dimensions::is_static), "Size must be known at compile time.");
	static std::size_t const value= dimensions::Num_Rows * dimensions::Num_Cols;
    };

    // return const-ref if matrix on stack and type itself if on heap 
    template <typename Matrix, bool on_stack>
    struct ref_on_stack
    {
	typedef Matrix type;
    };

    template <typename Matrix>
    struct ref_on_stack<Matrix, true>
    {
	typedef const Matrix& type;
    };

} // namespace detail

  
/// Dense matrix type
template <typename Value, typename Parameters = parameters<> >
class dense2D
    : public base_sub_matrix<Value, Parameters>,
    public mtl::detail::contiguous_memory_block< Value, Parameters::on_stack,
    detail::dense2D_array_size<Parameters, Parameters::on_stack>::value >,
    public crtp_base_matrix< dense2D<Value, Parameters>, Value, std::size_t >,
    public mat_expr< dense2D<Value, Parameters> >
{
    typedef dense2D                                           self;
    typedef base_sub_matrix<Value, Parameters>                super;
    typedef mtl::detail::contiguous_memory_block<Value, Parameters::on_stack,
	detail::dense2D_array_size<Parameters, Parameters::on_stack>::value>     memory_base;
    typedef mat_expr< dense2D<Value, Parameters> >            expr_base;
    typedef crtp_base_matrix< self, Value, std::size_t >      crtp_base;
    typedef crtp_matrix_assign< self, Value, std::size_t >    assign_base;
public:
    typedef Parameters                        parameters;
    typedef typename Parameters::orientation  orientation;
    typedef typename Parameters::index        index_type;
    typedef typename Parameters::dimensions   dim_type;
    typedef Value                             value_type;
    typedef const value_type&                 const_reference;
    typedef value_type&                       reference;

    typedef const value_type*                 const_pointer_type;
    typedef const_pointer_type                key_type;
    typedef std::size_t                       size_type;
    typedef dense_el_cursor<Value>            el_cursor_type;
    typedef dense2D_indexer                   indexer_type;

    // Self-similar type unless dimension is fixed
    // Not supported for the moment
    typedef self                              sub_matrix_type;

protected:
    // Obviously, the next 3 functions must be called after setting dimensions
    void set_nnz() { this->my_nnz = this->num_rows() * this->num_cols(); }
    void set_ldim(row_major) { ldim = this->num_cols(); }
    void set_ldim(col_major) { ldim = this->num_rows(); }
    void set_ldim() { set_ldim(orientation()); }

    void init()
    {
	set_nnz(); set_ldim(); // set_to_zero(*this);
    }

public:
    /// Default constructor, if compile time matrix size allocate memory
    dense2D() : memory_base(dim_type().num_rows() * dim_type().num_cols())
    {
	init();
    }

    /// Constructor that only sets dimensions, only for run-time dimensions
    explicit dense2D(mtl::non_fixed::dimensions d)
	: super(d), memory_base(d.num_rows() * d.num_cols())
    {
	init();
    }

    /// Most common constructor from number of rows and columns
    explicit dense2D(size_type num_rows, size_type num_cols)
	: super(dim_type(num_rows, num_cols)),
	memory_base(num_rows * num_cols)
    {
	init();
    }

    /// Constructor that sets dimensions and pointer to external data
    explicit dense2D(mtl::non_fixed::dimensions d, value_type* a)
	: super(d), memory_base(a, d.num_rows() * d.num_cols())
    {
	init();
    }

    /// Constructor that sets dimensions and pointer to external data
    explicit dense2D(size_type num_rows, size_type num_cols, value_type* a)
	: super(mtl::non_fixed::dimensions(num_rows, num_cols)), memory_base(a, num_rows * num_cols)
    {
	init();
    }

    /// Constructor for compile time matrix size
    /** sets dimensions and pointer to external data **/
    explicit dense2D(value_type* a)
	: super(), memory_base(a, dim_type().num_rows() * dim_type().num_cols())
    {
	MTL_STATIC_ASSERT((dim_type::is_static), "Size must be known at compile time.");
	init();
    }

    /// Default copy constructor
    dense2D(const self& m)
	: super(dim_type(m.num_rows(), m.num_cols())),
	memory_base(m)
    {
	// In case of sub-matrices we need m's ldim -> init doesn't work
	this->my_nnz = m.my_nnz; ldim = m.ldim;
    }

    /// Clone constructor, copies every source including sub-matrices and other matrices with references
    explicit dense2D(const self& m, clone_ctor)
	: super(mtl::non_fixed::dimensions(m.num_rows(), m.num_cols())),
	memory_base(m, clone_ctor())
    {
	init();
	*this = m;
    }

    /// General copy constructor, uses functionality from CRTP base
    template <typename MatrixSrc>
    explicit dense2D(const MatrixSrc& src)
	: super(), memory_base(dim_type().num_rows() * dim_type().num_cols())
    {
	init();
	*this = src;
    }

#ifdef MTL_WITH_MOVE
    /// Move constructor
    dense2D(self&& src)
	: super(std::move(src)), memory_base(std::move(src)), ldim(src.ldim)
    {}
#endif	

    /// Constructor for creating sub-matrices
    template <typename MatrixSrc>
    dense2D(MatrixSrc& matrix, dense2D_sub_ctor, 
	    size_type begin_r, size_type end_r, size_type begin_c, size_type end_c)
      : super(mtl::non_fixed::dimensions(matrix.num_rows(), matrix.num_cols())),
        memory_base(matrix.data, (end_r - begin_r) * (end_c - begin_c), true)
    {
	sub_matrix_constructor(matrix, begin_r, end_r, begin_c, end_c, boost::mpl::bool_<memory_base::on_stack>());
    }

#if defined(MTL_WITH_INITLIST) && defined(MTL_WITH_AUTO) && defined(MTL_WITH_RANGEDFOR)
    /// Constructor for initializer list \p values 
    template <typename Value2>
    dense2D(std::initializer_list<std::initializer_list<Value2> > values)
      : super(mtl::non_fixed::dimensions(values.size(), values.size()? values.begin()->size() : 0)),
    	memory_base(this->num_rows() * this->num_cols()) 
    {
    	init();
	*this= values;
    }
#endif

  private:
    template <typename MatrixSrc>
    void sub_matrix_constructor(MatrixSrc& matrix, size_type begin_r, size_type end_r, 
				size_type begin_c, size_type end_c, boost::mpl::false_)
    {
	matrix.check_ranges(begin_r, end_r, begin_c, end_c);
	
	if(end_r <= begin_r || end_c <= begin_c)
	    set_ranges(0, 0);
	else {
	    // Leading dimension doesn't change
	    this->data += matrix.indexer(matrix, begin_r, begin_c);  // Takes care of indexing
	    set_ranges(end_r - begin_r, end_c - begin_c);
	}
	this->my_nnz= matrix.nnz(); ldim= matrix.get_ldim();
    }

    template <typename MatrixSrc>
    void sub_matrix_constructor(MatrixSrc&, size_type, size_type, 
				size_type, size_type, boost::mpl::true_)
    {
	MTL_THROW(logic_error("Matrices cannot be used as sub-matrices!"));
    }

  public:

#ifdef MTL_WITH_MOVE
    /// Move assignment for data on heap
    self& operator=(self&& src)
    {	return self_assign(src, boost::mpl::bool_<memory_base::on_stack>());    }

    /// (Copy) Assignment
    self& operator=(const self& src)
    {	return self_assign(src, boost::mpl::true_());    }
    
#else   
    /// (Copy) Assignment
    self& operator=(typename detail::ref_on_stack<self, memory_base::on_stack>::type src)
    {
	return self_assign(src, boost::mpl::bool_<memory_base::on_stack>());
    }
#endif

  private:
    // Already copied for lvalues -> data can be stolen (need non-const ref)  
    self& self_assign(self& src, boost::mpl::false_)
    {
	// Self-copy would be an indication of an error
	assert(this != &src);
	// std::cout << "In move assignment: this* = \n" << *this << "src = \n" << src;

	this->checked_change_dim(src.num_rows(), src.num_cols());
	if (this->category == memory_base::view || src.category == memory_base::view)
	    matrix_copy(src, *this);
	else {
	    if (this->num_rows() != src.num_rows() || this->num_cols() != src.num_cols()) {
		super::change_dim(src.num_rows(), src.num_cols());
		init();
	    }
	    memory_base::move_assignment(src);
	}
	// std::cout << "End of move assignment: this* = \n" << *this;
	return *this;
    }

    // For matrices with data on stack (or lvalues in C++11)
    self& self_assign(const self& src, boost::mpl::true_)
    {
	if (this != &src) {
	    this->checked_change_dim(src.num_rows(), src.num_cols());
	    matrix_copy(src, *this);
	}
	return *this;
    }
  public:

    // import operators from CRTP base class
#if 0 // def __PGI
    using crtp_base::operator=;
#else
    using assign_base::operator=;
#endif

    /// Change dimension, can keep old data
    void change_dim(size_type r, size_type c, bool keep_data = false)
    {
	change_dim(r, c, keep_data, mtl::traits::is_static<self>());
    }

  private:
    void change_dim(size_type r, size_type c, bool keep_data, boost::mpl::false_)
    {
	if (r == this->num_rows() && c == this->num_cols())
	    return;

	self temp;
	if (keep_data) {
	    temp.super::change_dim(this->num_rows(), this->num_cols());
	    temp.init();
	    temp.memory_base::move_assignment(*this);
	}
	memory_base::realloc(r*c);
	super::change_dim(r, c);
	init();
	if (keep_data) {
	    if (r > temp.num_rows() || c > temp.num_cols()){
		set_to_zero(*this);
#if 0
		irange rr(0, std::min(r,temp.num_rows())), cr(0, std::min(c,temp.num_cols()));
		*this[rr][cr]= temp[rr][cr];
#endif
		sub_matrix(*this,0,std::min(r,temp.num_rows()),0,std::min(c,temp.num_cols()))
		    = sub_matrix(temp,0,std::min(r,temp.num_rows()),0,std::min(c,temp.num_cols()));
	    } else 
		*this = temp[irange(0, r)][irange(0, c)];
	}
    }

    void change_dim(size_type MTL_DEBUG_ARG(r), size_type MTL_DEBUG_ARG(c), bool, boost::mpl::true_)
    {	assert(r == this->num_rows() && c == this->num_cols());    }

 public:
    /// Check whether indices r and c are in range
    bool check_indices(size_t r, size_t c) const
    {	return r >= this->begin_row() && r < this->end_row() && c >= this->begin_col() && c < this->end_col();    }

    /// Constant access to element
    const_reference operator() (size_t r, size_t c) const 
    {
	MTL_CRASH_IF(is_negative(r) || r >= this->num_rows() || is_negative(c) || c >= this->num_cols(), 
		     "Index out of range!");
        return this->data[indexer(*this, r, c)];
    }

    /// Mutable access to element
    value_type& operator() (size_t r, size_t c)
    {
	MTL_CRASH_IF(is_negative(r) || r >= this->num_rows() || is_negative(c) || c >= this->num_cols(), 
		     "Index out of range!");
	return this->data[indexer(*this, r, c)]; 
    }    

    // offset regarding c-style indices
    size_t c_offset(size_t r, size_t c) const
    {	return indexer.offset(ldim, r, c, orientation());    }

    /// Get lower dimension [advanced]
    size_type get_ldim() const
    {	return ldim;    }

    /// Swap two matrices
    friend void swap(self& matrix1, self& matrix2)
    {
	swap(static_cast<memory_base&>(matrix1), static_cast<memory_base&>(matrix2));
	swap(static_cast<super&>(matrix1), static_cast<super&>(matrix2));
	std::swap(matrix1.ldim, matrix2.ldim);
    }

    void crop() {} ///< Delete structural zeros; only dummy here

    /// Address of first data entry (mutable); to be used with care. [advanced]
    value_type* address_data() { return this->data; }
    /// Address of first data entry (constant); to be used with care. [advanced]
    const value_type* address_data() const { return this->data; }

    /// Whether data is stored in strides
    bool has_strided_data() const { return this->category != this->own; }
    
  protected:
    
    // Set ranges from begin_r to end_r and begin_c to end_c
    void set_ranges(size_type begin_r, size_type end_r, size_type begin_c, size_type end_c)
    {
	super::set_ranges(begin_r, end_r, begin_c, end_c);
	set_nnz();
    }
	
    // Set ranges to a num_row x num_col matrix, keeps indexing
    void set_ranges(size_type num_rows, size_type num_cols)
    {
	set_ranges(this->begin_row(), this->begin_row() + num_rows, 
		   this->begin_col(), this->begin_col() + num_cols);
    }
    
  public:

    indexer_type  indexer;

    friend class dense2D_indexer;

#if !defined(_MSC_VER) || _MSC_VER != 1400 // Bug in MSVC 2005
    template <typename> friend struct sub_matrix_t;
    template <typename, typename> friend struct mtl::traits::range_generator;
    template <typename, typename, bool> friend struct mtl::traits::detail::dense2D_iterator_range_generator;

  protected:
#endif

    // Leading dimension is minor dimension in original matrix 
    // Opposed to other dims doesn't change in sub-matrices
    size_type     ldim; 

}; // dense2D


// ================
// Free functions
// ================


/// Number of rows
template <typename Value, typename Parameters>
typename dense2D<Value, Parameters>::size_type
inline num_rows(const dense2D<Value, Parameters>& matrix)
{
    return matrix.num_rows();
}

/// Number of columns
template <typename Value, typename Parameters>
typename dense2D<Value, Parameters>::size_type
inline num_cols(const dense2D<Value, Parameters>& matrix)
{
    return matrix.num_cols();
}

/// Size of the matrix, i.e. the number of row times columns
template <typename Value, typename Parameters>
typename dense2D<Value, Parameters>::size_type
inline size(const dense2D<Value, Parameters>& matrix)
{
    return matrix.num_cols() * matrix.num_rows();
}


}

using mat::dense2D;

} // namespace mtl::matrix


namespace mtl { namespace traits {


    // VC 8.0 finds ambiguity with mtl::tag::dense2D (I wonder why, especially here)
    using mtl::mat::dense2D;

    // ================
    // Range generators
    // For cursors
    // ================

    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::all, dense2D<Value, Parameters> >
      : detail::dense_element_range_generator<dense2D<Value, Parameters>,
					      dense_el_cursor<Value>, complexity_classes::linear_cached>
    {};

    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::nz, dense2D<Value, Parameters> >
      : detail::dense_element_range_generator<dense2D<Value, Parameters>,
					      dense_el_cursor<Value>, complexity_classes::linear_cached>
    {};

    namespace detail 
    {
	// complexity of dense row cursor depends on storage scheme
	// if orientation is row_major then complexity is cached_linear, otherwise linear
	template <typename Orientation> struct dense2D_rc {};
	template<> struct dense2D_rc<row_major>
	{
	    typedef complexity_classes::linear_cached type;
	};
	template<> struct dense2D_rc<col_major>
	{
	    typedef complexity_classes::linear type;
	};

	// Complexity of column cursor is of course opposite
	template <typename Orientation> struct dense2D_cc
	    : dense2D_rc<typename mtl::traits::transposed_orientation<Orientation>::type>
	{};
    }

    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::row, dense2D<Value, Parameters> >
	: detail::all_rows_range_generator<dense2D<Value, Parameters>, 
					   typename detail::dense2D_rc<typename Parameters::orientation>::type>
    {};
 
    // For a cursor pointing to some row give the range of elements in this row 
    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::nz, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, glas::tag::row, 2> >
    {
	typedef dense2D<Value, Parameters>                                            matrix;
	typedef typename matrix::size_type                                            size_type;
	typedef detail::sub_matrix_cursor<matrix, glas::tag::row, 2>               cursor;

	// linear for col_major and linear_cached for row_major
	typedef typename detail::dense2D_rc<typename Parameters::orientation>::type   complexity;
	static int const                                                              level = 1;

	typedef typename boost::mpl::if_<
	    boost::is_same<typename Parameters::orientation, row_major>
	  , dense_el_cursor<Value>
	  , strided_dense_el_cursor<Value>
	>::type type;  

      private:

	type dispatch(cursor const& c, size_type col, row_major) const
	{
	    return type(c.ref, c.key, col);
	}
	type dispatch(cursor const& c, size_type col, col_major) const
	{
	    return type(c.ref, c.key, col, c.ref.ldim);
	}

      public:

	type begin(cursor const& c) const
	{
	    return dispatch(c, c.ref.begin_col(), typename matrix::orientation());
	}
	type end(cursor const& c) const
	{
	    return dispatch(c, c.ref.end_col(), typename matrix::orientation());
	}	
	type lower_bound(cursor const& c, size_type position) const
	{
	    return dispatch(c, std::min(c.ref.end_col(), position), typename matrix::orientation());
	}
    };

    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::all, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, glas::tag::row, 2> >
        : range_generator<glas::tag::nz, 
			  detail::sub_matrix_cursor<dense2D<Value, Parameters>, glas::tag::row, 2> >
    {};


    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::col, dense2D<Value, Parameters> >
	: detail::all_cols_range_generator<dense2D<Value, Parameters>, 
					   typename detail::dense2D_cc<typename Parameters::orientation>::type>
    {};
 
    // For a cursor pointing to some row give the range of elements in this row 
    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::nz, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, glas::tag::col, 2> >
    {
	typedef dense2D<Value, Parameters>                                            matrix;
	typedef typename matrix::size_type                                            size_type;
	typedef detail::sub_matrix_cursor<matrix, glas::tag::col, 2>               cursor;	
	typedef typename detail::dense2D_cc<typename Parameters::orientation>::type   complexity;
	static int const                                                              level = 1;

	typedef typename boost::mpl::if_<
	    boost::is_same<typename Parameters::orientation, col_major>
	  , dense_el_cursor<Value>
	  , strided_dense_el_cursor<Value>
	>::type type;  

      private:
	type dispatch(cursor const& c, size_type row, col_major) const
	{
	    return type(c.ref, row, c.key);
	}
	type dispatch(cursor const& c, size_type row, row_major) const
	{
	    return type(c.ref, row, c.key, c.ref.ldim);
	}

      public:
	type begin(cursor const& c) const
	{
	    return dispatch(c, c.ref.begin_row(), typename matrix::orientation());
	}
	type end(cursor const& c) const
	{
	    return dispatch(c, c.ref.end_row(), typename matrix::orientation());
	}
	type lower_bound(cursor const& c, size_type position) const
	{
	    return dispatch(c, std::min(c.ref.end_row(), position), typename matrix::orientation());
	}	
    };

    template <typename Value, typename Parameters>
    struct range_generator<glas::tag::all, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, glas::tag::col, 2> >
      : public range_generator<glas::tag::nz, 
			       detail::sub_matrix_cursor<dense2D<Value, Parameters>, glas::tag::col, 2> >
    {};

// =============
// For iterators
// =============


    namespace detail {

        // Traversal along major dimension first and then along minor
        template <typename OuterTag, typename Orientation>
        struct major_traversal
        {
	    static const bool value= false;
        };
          
        template <> struct major_traversal<glas::tag::row, row_major>
        {
	    static const bool value= true;
        };
        
        template <> struct major_traversal<glas::tag::col, col_major>
        {
	    static const bool value= true;
        };


        template <typename OuterTag, typename Matrix, bool is_const>
        struct dense2D_iterator_range_generator
        {
	    typedef Matrix                                                                matrix_type;
	    typedef typename matrix_type::size_type                                       size_type;
	    typedef typename matrix_type::value_type                                      value_type;
	    typedef typename matrix_type::parameters                                      parameters;
	    typedef detail::sub_matrix_cursor<matrix_type, OuterTag, 2>                   cursor;

	    // if traverse first along major dimension then memory access is contiguous (otherwise strided)
	    typedef typename boost::mpl::if_<
		major_traversal<OuterTag, typename parameters::orientation> 
	      , complexity_classes::linear_cached
	      , complexity_classes::linear
	    >::type                                                                       complexity;
	    static int const                                                              level = 1;

	    // if traverse first along major dimension use pointer otherwise strided iterator
	    typedef typename boost::mpl::if_<
		major_traversal<OuterTag, typename parameters::orientation> 
	      , typename boost::mpl::if_c<
    	            is_const 
		  , const value_type*
    	          , value_type*
		>::type
	      , typename boost::mpl::if_c<
    	            is_const 
		  , strided_dense_el_const_iterator<value_type>
    	          , strided_dense_el_iterator<value_type>
    	        >::type
	    >::type type;  

        private:
	    // if traverse first along major dim. then return address as pointer
	    type dispatch(cursor const& c, size_type row, size_type col, complexity_classes::linear_cached) const
	    {
		matrix_type& ma= const_cast<matrix_type&>(c.ref);
		return ma.elements() + ma.indexer(ma, row, col); // &ref[row][col];
	    }

	    // otherwise strided 
	    type dispatch(cursor const& c, size_type row, size_type col, complexity_classes::linear) const
	    {
		// cast const away (is dirty and should be improved later (cursors must distinct constness))
		matrix_type& ref= const_cast<matrix_type&>(c.ref);
		return type(ref, row, col, ref.ldim);
	    }

	    type begin_dispatch(cursor const& c, glas::tag::row) const
	    {
		return dispatch(c, c.key, c.ref.begin_col(), complexity());
	    }
	    
	    type end_dispatch(cursor const& c, glas::tag::row) const
	    {
		return dispatch(c, c.key, c.ref.end_col(), complexity());
	    }


	    type begin_dispatch(cursor const& c, glas::tag::col) const
	    {
		return dispatch(c, c.ref.begin_row(), c.key, complexity());
	    }

	    type end_dispatch(cursor const& c, glas::tag::col) const
	    {
		return dispatch(c, c.ref.end_row(), c.key, complexity());
	    }

        public:

	    type begin(cursor const& c) const
	    {
		return begin_dispatch(c, OuterTag());
	    }

	    type end(cursor const& c) const
	    {
		return end_dispatch(c, OuterTag());
	    }	
        };

    } // namespace detail

        
    template <typename Value, typename Parameters, typename OuterTag>
    struct range_generator<tag::iter::nz, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, OuterTag, 2> >
      : public detail::dense2D_iterator_range_generator<OuterTag, dense2D<Value, Parameters>, false>
    {};

    template <typename Value, typename Parameters, typename OuterTag>
    struct range_generator<tag::iter::all, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, OuterTag, 2> >
      : public detail::dense2D_iterator_range_generator<OuterTag, dense2D<Value, Parameters>, false>
    {};

    template <typename Value, typename Parameters, typename OuterTag>
    struct range_generator<tag::const_iter::nz, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, OuterTag, 2> >
      : public detail::dense2D_iterator_range_generator<OuterTag, dense2D<Value, Parameters>, true>
    {};

    template <typename Value, typename Parameters, typename OuterTag>
    struct range_generator<tag::const_iter::all, 
			   detail::sub_matrix_cursor<dense2D<Value, Parameters>, OuterTag, 2> >
      : public detail::dense2D_iterator_range_generator<OuterTag, dense2D<Value, Parameters>, true>
    {};


}} // namespace mtl::traits

namespace mtl { namespace mat {

    // ==========
    // Sub matrix
    // ==========

    template <typename Value, typename Parameters>
    struct sub_matrix_t<dense2D<Value, Parameters> >
    {
        typedef dense2D<Value, Parameters>      matrix_type;
	// copy orientation, ignore index, set dimension to non-fixed and on_stack to false
	typedef parameters<typename Parameters::orientation> para_type; 

        typedef dense2D<Value, para_type>       sub_matrix_type;
        typedef sub_matrix_type const           const_sub_matrix_type;
        typedef typename matrix_type::size_type size_type;
        
        sub_matrix_type operator()(matrix_type& matrix, size_type begin_r, size_type end_r, size_type begin_c, size_type end_c)
        {
	    return sub_matrix_type(matrix, dense2D_sub_ctor(), begin_r, end_r, begin_c, end_c);
        }

        const_sub_matrix_type
        operator()(matrix_type const& matrix, size_type begin_r, size_type end_r, size_type begin_c, size_type end_c)
        {
	    // To minimize code duplication, we use the non-const version
	    sub_matrix_type tmp((*this)(const_cast<matrix_type&>(matrix), begin_r, end_r, begin_c, end_c));
	    return tmp;
        }	
    };
        
}} // mtl::matrix

namespace mtl {

    // Enable cloning of dense matrices
    template <typename Value, typename Parameters>
    struct is_clonable< mtl::mat::dense2D<Value, Parameters> > : boost::mpl::true_ {};
        
} // namespace mtl



namespace math {

    // Multiplicative identities of matrices
    template <typename Value, typename Parameters>
    struct identity_t< mult<mtl::mat::dense2D<Value, Parameters> >, mtl::mat::dense2D<Value, Parameters> >
        : public binary_function< mult<mtl::mat::dense2D<Value, Parameters> >, 
				       mtl::mat::dense2D<Value, Parameters>, 
				       mtl::mat::dense2D<Value, Parameters> >
    {
        typedef mtl::mat::dense2D<Value, Parameters>  matrix_type;

        matrix_type operator() (const mult<matrix_type>&, const matrix_type& ref) const
        {
	    matrix_type tmp(ref);
	    tmp= one(typename matrix_type::value_type());
	    return tmp;
        }
    };
        
} // namespace math


#endif // MTL_DENSE2D_INCLUDE


/*
Limitations:
- with compile-time constant dimension, submatrices are not supported (would violate self-similarity)
- Element cursor doesn't work for sub-matrices (not contiguous)
*/
