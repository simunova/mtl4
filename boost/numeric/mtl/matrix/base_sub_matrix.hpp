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

#ifndef MTL_BASE_SUB_MATRIX_INCLUDE
#define MTL_BASE_SUB_MATRIX_INCLUDE

#include <algorithm>
#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/detail/index.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>

namespace mtl { namespace mat {

// Base class for sub-matrices
// Contains only very simple functionality that is used in all sub-matrices
// But also used in some complete matrices
template <class Elt, class Parameters>
struct base_sub_matrix 
{
    typedef Elt                               value_type;
    typedef typename Parameters::orientation  orientation;
    typedef typename Parameters::index        index_type;
    typedef typename Parameters::dimensions   dim_type;
    static bool const                         on_stack= Parameters::on_stack;
    typedef std::size_t                       size_type;
    typedef base_sub_matrix                   self;

  protected:
    size_type                       my_nnz,       // # of non-zeros, to be set by derived matrix (drop maybe?)
                                    my_begin_row, my_end_row,
                                    my_begin_col, my_end_col;

    void constructor_helper(const dim_type& dim)
    {
	my_begin_row= index::change_to(index_type(), 0);
	my_end_row=   index::change_to(index_type(), dim.num_rows());
	my_begin_col= index::change_to(index_type(), 0);
	my_end_col=   index::change_to(index_type(), dim.num_cols());
	my_nnz= 0;
    }

  public:
    // base_sub_matrix() :  my_nnz(0), my_begin_row(0), my_end_row(0), my_begin_col(0), my_end_col(0) {}
   
    base_sub_matrix() 
    {
	// With no static dimension information, it is by default 0
	constructor_helper(dim_type());
    }

    explicit base_sub_matrix(const dim_type& d) 
    //explicit base_sub_matrix(mtl::non_fixed::dimensions d) 
    {
	constructor_helper(d);
    }

    friend void swap(self& x, self& y)
    {
	std::swap(x.my_nnz, y.my_nnz);
	std::swap(x.my_begin_row, y.my_begin_row);
	std::swap(x.my_end_row, y.my_end_row);
	std::swap(x.my_begin_col, y.my_begin_col);
	std::swap(x.my_end_col, y.my_end_col);
    }

    // Either changed matrix is uninitialized (i.e. 0x0) or dimensions are equal
    void check_dim(size_type MTL_DEBUG_ARG(num_rows), size_type MTL_DEBUG_ARG(num_cols) ) const
    {
	MTL_CRASH_IF(this->num_rows() * this->num_cols() != 0
		   && (this->num_rows() != num_rows || this->num_cols() != num_cols),
		   "Incompatible size");
    }

protected:
    void change_dim(non_fixed::dimensions d) { constructor_helper(d); }    
    template <std::size_t Rows, std::size_t Cols>
    void change_dim(fixed::dimensions<Rows, Cols> d) { check_dim(d.num_rows(), d.num_cols()); } 


   
    void change_dim(size_type r, size_type c) {	change_dim(dim_type(r, c));  }    

    void set_ranges(size_type br, size_type er, size_type bc, size_type ec)
    {
	MTL_CRASH_IF(br > er, "begin row > end row");
	MTL_CRASH_IF(bc > ec, "begin column > end column");
	my_begin_row= br; my_end_row= er; my_begin_col= bc; my_end_col= ec;
    }

public:
    void check_ranges(size_type MTL_DEBUG_ARG(begin_r), size_type MTL_DEBUG_ARG(end_r), 
		      size_type MTL_DEBUG_ARG(begin_c), size_type MTL_DEBUG_ARG(end_c) ) const
    {
	MTL_CRASH_IF(begin_r < begin_row(), "begin_row out of range");
	// if (end_r > end_row()) std::cout << "end_row out of range\n";
	MTL_CRASH_IF(end_r > end_row(), "end_row out of range");
			      
	MTL_CRASH_IF(begin_c < begin_col(), "begin_col out of range");
	MTL_CRASH_IF(end_c > end_col(), "end_col out of range");
    }

    explicit base_sub_matrix(size_type br, size_type er, size_type bc, size_type ec) : my_nnz(0)
    {
	set_ranges(br, er, bc, ec);
    }
 
    // Number of rows
    size_type num_rows() const 
    {
	return my_end_row - my_begin_row;
    }

    // First row taking indexing into account (already stored as such)
    size_type begin_row() const 
    {
	return my_begin_row;
    }
    
    // Past-end row taking indexing into account (already stored as such)
    size_type end_row() const 
    {
	return my_end_row;
    }

    // Number of columns
    size_type num_cols() const 
    {
	return my_end_col - my_begin_col;
    }

    // First column taking indexing into account (already stored as such)
    size_type begin_col() const 
    {
	return my_begin_col;
    }
    
    // Past-end column taking indexing into account (already stored as such)
    size_type end_col() const 
    {
	return my_end_col;
    }

    // Number of non-zeros
    size_type nnz() const
    {
      return my_nnz;
    }

  protected:
    // dispatched functions for major dimension
    size_type dim1(row_major) const 
    {
	return num_rows();
    }

    size_type dim1(col_major) const 
    {
	return num_cols();
    }

    // dispatched functions for minor dimension
    size_type dim2(row_major) const 
    {
	return num_cols();
    }

    size_type dim2(col_major) const 
    {
	return num_rows();
    }
  
    // Dispatched functions for major
    // Trailing _ due to conflicts with macro major
    size_type major_(size_type r, size_type, row_major) const
    {
	return r; 
    }

    size_type major_(size_type, size_type c, col_major) const
    {
	return c; 
    }    

  public:
    // return major dimension
    size_type dim1() const 
    {
	return dim1(orientation());
    }

    // return minor dimension
    size_type dim2() const 
    {
	return dim2(orientation());
    }

    // Returns the row for row_major otherwise the column
    // Trailing _ due to conflicts with macro major
    size_type major_(size_type r, size_type c) const
    {
	return major_(r, c, orientation());
    }

    // Returns the row for col_major otherwise the column
    // Trailing _ for consistency with major
    size_type minor_(size_type r, size_type c) const
    {
	return major_(c, r, orientation());
    }

    dim_type get_dimensions() const 
    {
	return dim_type(num_rows(), num_cols()); 
    }    

};


}} // namespace mtl::matrix

#endif // MTL_BASE_SUB_MATRIX_INCLUDE


/* 
   Question:
   - Shall we keep the position in the original matrix?
*/
