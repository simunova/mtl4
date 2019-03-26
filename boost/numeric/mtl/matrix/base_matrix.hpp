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

#ifndef MTL_BASE_MATRIX_INCLUDE
#define MTL_BASE_MATRIX_INCLUDE

#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/detail/index.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>

namespace mtl { namespace mat {
  
/// Base class for other matrices, contains only very simple functionality that is used in all matrices.
template <class Elt, class Parameters>
struct base_matrix 
{
    typedef base_matrix                       self;
    typedef Elt                               value_type;
    typedef typename Parameters::orientation  orientation;
    typedef typename Parameters::index        index_type;
    typedef typename Parameters::dimensions   dim_type;
    static bool const                         on_stack= Parameters::on_stack;
    typedef typename Parameters::size_type    size_type;
  protected:
    dim_type                        dim;       ///< # of rows and columns
    size_type                       my_nnz;    ///< # of non-zeros, to be set by derived matrix
    typedef mtl::traits::is_static<dim_type>       static_bool;
    
  public:
    base_matrix(size_type n= 0) :  my_nnz(n) {}

    /// Setting dimension
    explicit base_matrix(mtl::non_fixed::dimensions d, size_type n= 0) : dim(d), my_nnz(n) {}

	///	Swap base matrix
    friend void swap(self& x, self& y)
    {
	using std::swap;
	swap(x.my_nnz, y.my_nnz);
	swap(x.dim, y.dim);
    }

    /// Either matrix to be changed is uninitialized (i.e. 0x0) or dimensions are equal
    /** The matrices with dimension 0 x 0 are considered like stem cells: they can still
	change into an arbitrary dimension and are compatible with any other matrix.  Once a matrix has a non-trivial dimension
	it can be only changed explicitly and is only compatible with matrices of the same dimensionality. **/
    void check_dim(size_type MTL_DEBUG_ARG(num_rows), size_type MTL_DEBUG_ARG(num_cols)) const
    {
	MTL_CRASH_IF(this->num_rows() * this->num_cols() != 0
		   && (this->num_rows() != num_rows || this->num_cols() != num_cols),
		   "Incompatible size");
    }

#if 0
    /** Will fail for fixed::dimension **/
    void change_dim(mtl::non_fixed::dimensions d) { dim= d; }

    template <std::size_t Rows, std::size_t Cols> 
    void change_dim(mtl::fixed::dimensions<Rows, Cols> d) {}
#endif

    void change_dim(size_type r, size_type c, boost::mpl::false_) { dim= dim_type(r, c); }    
    void change_dim(size_type r, size_type c, boost::mpl::true_) { check_dim(r, c); }    

    void change_dim(size_type r, size_type c) {	change_dim(r, c, static_bool()); }    

public:
    /// Number of rows
    size_type num_rows() const 
    {
      return size_type(dim.num_rows()); // prevent narrowing warnings
    }
    /// First row taking indexing into account
    size_type begin_row() const 
    {
      return index::change_to(index_type(), 0);
    }
    /// Past-end row taking indexing into account
    size_type end_row() const 
    {
      return index::change_to(index_type(), num_rows());
    }

    /// number of colums
    size_type num_cols() const 
    {
		return size_type(dim.num_cols()); // prevent narrowing warnings
    }
    /// First column taking indexing into account
    size_type begin_col() const 
    {
      return index::change_to(index_type(), 0);
    }
    /// Past-end column taking indexing into account
    size_type end_col() const 
    {
      return index::change_to(index_type(), num_cols());
    }

    /// Number of non-zeros
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
    /// Major dimension
    size_type dim1() const 
    {
      return dim1(orientation());
    }

    /// Minor dimension
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
	
    // returns copy of dim
    dim_type get_dimensions() const 
    {
      return dim; 
    }    
};



}} // namespace mtl::matrix

#endif // MTL_BASE_MATRIX_INCLUDE
