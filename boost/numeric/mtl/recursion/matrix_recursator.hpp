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

#ifndef MTL_MATRIX_RECURATOR_INCLUDE
#define MTL_MATRIX_RECURATOR_INCLUDE

#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/sub_matrix.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/recursion/dim_splitter.hpp>
#include <boost/numeric/mtl/recursion/utility.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl { namespace mat {


template <typename Recursator1, typename Recursator2>
void inline equalize_depth(Recursator1& r1, Recursator2& r2);

template <typename Recursator1, typename Recursator2, typename Recursator3>
void inline equalize_depth(Recursator1& r1, Recursator2& r2, Recursator3& r3);


/*! Class for matrix recursator

    How to use this class is described in the \ref rec_intro "recursion introduction".

    \sa \ref mtl::mat::north_west, \ref mtl::mat::north_east, 
    \ref mtl::mat::south_west, \ref mtl::mat::south_east, 
    \ref mtl::mat::is_empty(const recursator<Matrix>&), 
    \ref mtl::mat::is_full(const recursator<Matrix>&), 
    \ref mtl::mat::num_rows(const recursator<Matrix>&), 
    \ref mtl::mat::num_cols(const recursator<Matrix>&),
    \ref mtl::mat::size(const recursator<Matrix>&)
**/
template <typename Matrix>
struct recursator
{
    typedef recursator                                     self;
    typedef Matrix                                                matrix_type;
    typedef typename sub_matrix_t<Matrix>::sub_matrix_type        sub_matrix_type;
    typedef typename sub_matrix_t<Matrix>::const_sub_matrix_type  const_sub_matrix_type;
    typedef typename Collection<Matrix>::size_type                size_type;
    typedef typename Collection<Matrix>::value_type               matrix_value_type;
    typedef recursion::outer_bound_splitter<self>                 splitter_type;

private:
    
    template <typename MatrixType> // why was it templated ???
    sub_matrix_type constructor_helper(MatrixType const& matrix)
    {
	return sub_matrix(matrix, matrix.begin_row(), matrix.end_row(),
			  matrix.begin_col(), matrix.end_col());
    }

    // For views without own data, we need to generate a new sub_matrix as shared_ptr
    template <typename MatrixType>
    sub_matrix_type constructor_helper(transposed_view<MatrixType> const& view)
    {
	typedef typename boost::remove_const<MatrixType>::type   tmp_type;
	typedef typename sub_matrix_t<tmp_type>::sub_matrix_type ref_sub_type;
	typedef boost::shared_ptr<ref_sub_type>                  pointer_type;
	typedef typename transposed_view<MatrixType>::other      ref_type;

	// Submatrix of referred matrix, colums and rows interchanged
	// Create a submatrix, whos address will be kept by transposed_view
	pointer_type p(new ref_sub_type(sub_matrix(const_cast<ref_type&>(view.ref), view.begin_col(), view.end_col(), 
						   view.begin_row(), view.end_row())));
	return sub_matrix_type(p); 
    }

public:
    /*! Construct a recursator from a matrix.
        \param matrix The matrix to which the recursator refers.
	\param bound  Explicit bound declaration; must not be smaller than the numbers of rows and the number of columns;
	              must also be a power of 2.

        Constructor takes the entire matrix as sub-matrix.
        This allows to have different type for the matrix and the sub-matrix.
    **/
    explicit recursator(Matrix const& matrix,
			       size_type bound= 0
			       ) 
	: my_sub_matrix(constructor_helper(matrix)), my_bound(recursion::outer_bound(matrix)),
	  my_first_row(0), my_first_col(0)         // splitter(*this)
    {
      if (bound == 0)
	my_bound= recursion::outer_bound(matrix);
      else {
	MTL_DEBUG_THROW_IF(!recursion::is_power_of_2(bound), range_error("Bound must be a power of 2"));
	MTL_DEBUG_THROW_IF(bound < num_rows(matrix) || bound < num_cols(matrix), 
			   range_error("Bound must not be smaller than matrix dimensions"));
	my_bound= bound;
      }
    }


private:

    template <typename SubMatrix>
    sub_matrix_type get_value_dispatch(const SubMatrix& , 
				   size_type br, size_type er, size_type bc, size_type ec) const
    {
	return sub_matrix(my_sub_matrix, br, er, bc, ec);
    }

    template <typename SubMatrix>
    sub_matrix_type get_value_dispatch(transposed_view<SubMatrix> view, 
				       size_type br, size_type er, size_type bc, size_type ec) const
    {
	typedef typename sub_matrix_t<SubMatrix>::sub_matrix_type   ref_sub_type;
	typedef boost::shared_ptr<ref_sub_type>                     pointer_type;
	typedef typename transposed_view<SubMatrix>::other          ref_type;

	pointer_type p(new ref_sub_type(sub_matrix(const_cast<ref_type&>(view.ref), bc, ec, br, er)));
	return sub_matrix_type(p); 	
    }


public:
    sub_matrix_type get_value() const
    {
	using std::min;
	size_type begin_row= my_sub_matrix.begin_row() + my_first_row,
	          end_row= min(begin_row + my_bound, my_sub_matrix.end_row()),
	          begin_col= my_sub_matrix.begin_col() + my_first_col,
	          end_col= min(begin_col + my_bound, my_sub_matrix.end_col());

#if 0
	std::cout << "get_value [" << begin_row << "-" << end_row << "]["
		  << begin_col << "-" << end_col << "]\n";
#endif
	return get_value_dispatch(my_sub_matrix, begin_row, end_row, begin_col, end_col);
    }

    /// Compute the sub-matrix corresponding to this recursator.
    sub_matrix_type operator*() const
    {
	return get_value();
    }

    // Returning quadrants for non-const recursator

    self north_west() const
    {
	self tmp(*this);
	tmp.my_bound >>= 1; // divide by 2
	return tmp;
    }

    self south_west() const
    {
	self tmp(*this);
	tmp.my_bound >>= 1; // divide by 2
	tmp.my_first_row += tmp.my_bound;
	return tmp;
    }

    self north_east() const
    {
	self tmp(*this);
	tmp.my_bound >>= 1; // divide by 2
	tmp.my_first_col += tmp.my_bound;
	return tmp;
    }

    self south_east() const
    {
	self tmp(*this);
	tmp.my_bound >>= 1; // divide by 20
	tmp.my_first_row += tmp.my_bound;
	tmp.my_first_col += tmp.my_bound;
	return tmp;
    }

    bool is_empty() const
    {
	return my_first_row >= num_rows(my_sub_matrix) || my_first_col >= num_cols(my_sub_matrix);
    }


    /// Return the bound of the recursator
    size_type bound() const
    {
	return my_bound;
    }

    /*! Set the bound of the recursator.
	\param b  The new virtual bound; must be a power of 2.

        This function allows to declare a virtual bound smaller than the number of rows and/or columns.
	It must be used with uttermost care.
    **/
    void set_bound(size_type b)
    {
	my_bound= b;
    }

    template <typename R1, typename R2> friend void equalize_depth (R1&, R2&);   
    template <typename R1, typename R2, typename R3> friend void equalize_depth (R1&, R2&, R3&);

    template <typename M> friend typename recursator<M>::size_type num_rows(const recursator<M>& rec);
    template <typename M> friend typename recursator<M>::size_type num_cols(const recursator<M>& rec);

    // Dirty feature to be used with care
    matrix_value_type* first_address()
    {
	return &my_sub_matrix[my_first_row][my_first_col];
    }

    const matrix_value_type* first_address() const
    {
	return &my_sub_matrix[my_first_row][my_first_col];
    }
    
  protected:
    sub_matrix_type     my_sub_matrix; /// Referred matrix (from which the sub-matrices are built)
    size_type           my_bound,      /// Virtual matrix size, i.e. upper bound for size of sub-matrix.
	                my_first_row,  /// Row of first entry in submatrix
                        my_first_col;  /// Row of first entry in submatrix 
};

#if 0

// Obsolete, only left in code because discussed in a paper
// To use recursator with const matrices Reference must be 'Matrix const&'
template <typename Matrix, typename Splitter = recursion::max_dim_splitter<Matrix> >
struct recursator_s
{
    typedef recursator_s                                    self;
    typedef Matrix                                                matrix_type;
    typedef Splitter                                              splitter_type;
    typedef typename sub_matrix_t<Matrix>::sub_matrix_type        sub_matrix_type;
    typedef typename sub_matrix_t<Matrix>::const_sub_matrix_type  const_sub_matrix_type;
    typedef typename Matrix::size_type                            size_type;
    // typedef outer_bound_splitter<self>                            splitter_type;

private:
    
    // template <typename Matrix> why was it templated ???
    sub_matrix_type constructor_helper(Matrix const& matrix)
    {
	return sub_matrix(matrix, matrix.begin_row(), matrix.end_row(),
			  matrix.begin_col(), matrix.end_col());
    }

    // For views without own data, we need to generate a new sub_matrix as shared_ptr
    // template <typename Matrix>
    sub_matrix_type constructor_helper(transposed_view<Matrix> const& matrix)
    {
	typedef typename sub_matrix_t<Matrix>::sub_matrix_type   ref_sub_type;
	typedef boost::shared_ptr<ref_sub_type>                  pointer_type;

	// Submatrix of referred matrix, colums and rows interchanged
	// Create a submatrix, whos address will be kept by transposed_view
	pointer_type p(new ref_sub_type(sub_matrix(matrix.ref, matrix.begin_col(), matrix.end_col(), 
						   matrix.begin_row(), matrix.end_row())));
	return sub_matrix_type(p); 
    }

public:
    // Constructor takes the whole matrix as sub-matrix
    // This allows to have different type for the matrix and the sub-matrix
    // This also enables matrices to have references as sub-matrices
    explicit recursator_s(Matrix const& matrix, size_type bound= 0) 
	: my_sub_matrix(constructor_helper(matrix)), my_bound(outer_bound(matrix)),
	  splitter(my_sub_matrix)
    {
      if (bound == 0)
	my_bound= outer_bound(matrix);
      else {
	assert(is_power_of_2(bound));
	assert(bound >= matrix.num_rows() && bound >= matrix.num_cols());
	my_bound= bound;
      }
    }

    // Sub-matrices are copied directly
    // explicit recursator(sub_matrix_type sub_matrix) : my_sub_matrix(sub_matrix) {}
    
    sub_matrix_type& get_value()
    {
	return my_sub_matrix;
    }

    sub_matrix_type const& get_value() const
    {
	return my_sub_matrix;
    }

    // Returning quadrants for non-const recursator

    self north_west()
    {
	sub_matrix_type sm(sub_matrix(my_sub_matrix, my_sub_matrix.begin_row(), splitter.row_split(),
				      my_sub_matrix.begin_col(), splitter.col_split()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    self south_west()
    {
	sub_matrix_type sm(sub_matrix(my_sub_matrix, splitter.row_split(), my_sub_matrix.end_row(), 
				      my_sub_matrix.begin_col(), splitter.col_split()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    self north_east()
    {
	sub_matrix_type sm(sub_matrix(my_sub_matrix, my_sub_matrix.begin_row(), splitter.row_split(),
				      splitter.col_split(), my_sub_matrix.end_col()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    self south_east()
    {
	sub_matrix_type sm(sub_matrix(my_sub_matrix, splitter.row_split(), my_sub_matrix.end_row(), 
				      splitter.col_split(), my_sub_matrix.end_col()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    // Returning quadrants for const recursator

    self const north_west() const
    {
	sub_matrix_type sm(sub_matrix(const_cast<self*>(this)->my_sub_matrix, my_sub_matrix.begin_row(), splitter.row_split(),
				      my_sub_matrix.begin_col(), splitter.col_split()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    self const south_west() const 
    {
	sub_matrix_type sm(sub_matrix(const_cast<self*>(this)->my_sub_matrix, splitter.row_split(), my_sub_matrix.end_row(), 
				      my_sub_matrix.begin_col(), splitter.col_split()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    self const north_east() const 
    {
	sub_matrix_type sm(sub_matrix(const_cast<self*>(this)->my_sub_matrix, my_sub_matrix.begin_row(), splitter.row_split(),
				      splitter.col_split(), my_sub_matrix.end_col()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    self const south_east() const 
    {
	sub_matrix_type sm(sub_matrix(const_cast<self*>(this)->my_sub_matrix, splitter.row_split(), my_sub_matrix.end_row(), 
				      splitter.col_split(), my_sub_matrix.end_col()));
	self tmp(sm, my_bound / 2);
	return tmp;
    }

    // Checking whether a quadrant is empty

    // For completeness
    bool north_west_empty() const
    {
	return false;
    }

    bool north_east_empty() const
    {
	return splitter.col_split() == my_sub_matrix.end_col();
    }

    bool south_west_empty() const
    {
	return splitter.row_split() == my_sub_matrix.end_row();
    }

    bool south_east_empty() const
    {
	return splitter.row_split() == my_sub_matrix.end_row() 
	       || splitter.col_split() == my_sub_matrix.end_col();
    }

    bool is_empty() const
    {
	return my_sub_matrix.begin_row() == my_sub_matrix.end_row()
	       || my_sub_matrix.begin_col() == my_sub_matrix.end_col();
    }

#if 0
    bool is_leaf() const
    {
	return my_sub_matrix.num_rows() < 2 || my_sub_matrix.num_cols() < 2;
    }
#endif

    size_type bound() const
    {
	assert(my_bound >= my_sub_matrix.num_rows() && my_bound >= my_sub_matrix.num_cols());
	return my_bound;
    }

    template <typename R1, typename R2> friend void equalize_depth (R1&, R2&);   
    template <typename R1, typename R2, typename R3> friend void equalize_depth (R1&, R2&, R3&);

  protected:
    sub_matrix_type     my_sub_matrix;
    size_type           my_bound;
    splitter_type       splitter;
};

#endif

template <typename Recursator1, typename Recursator2>
void inline equalize_depth(Recursator1& r1, Recursator2& r2)
{
    typename Recursator1::size_type max_bound= std::max(r1.bound(), r2.bound());
    r1.my_bound= max_bound;
    r2.my_bound= max_bound;
}

template <typename Recursator1, typename Recursator2, typename Recursator3>
void inline equalize_depth(Recursator1& r1, Recursator2& r2, Recursator3& r3)
{
    typename Recursator1::size_type max_bound= std::max(std::max(r1.bound(), r2.bound()), r3.bound());
    r1.my_bound= max_bound;
    r2.my_bound= max_bound;
    r3.my_bound= max_bound;
}


// Define free functions (from member functions)

/*! Compute the north-west quadrant of a recursator (i.e. its referred matrix).
    The result is itself a recursator.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
recursator<Matrix> inline north_west(const recursator<Matrix>& rec)
{
    return rec.north_west();
}

/*! Compute the north-east quadrant of a recursator (i.e. its referred matrix).
    The result is itself a recursator.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
recursator<Matrix> inline north_east(const recursator<Matrix>& rec)
{
    return rec.north_east();
}

/*! Compute the south-west quadrant of a recursator (i.e. its referred matrix).
    The result is itself a recursator.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
recursator<Matrix> inline south_west(const recursator<Matrix>& rec)
{
    return rec.south_west();
}

/*! Compute the south-east quadrant of a recursator (i.e. its referred matrix).
    The result is itself a recursator.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
recursator<Matrix> inline south_east(const recursator<Matrix>& rec)
{
    return rec.south_east();
}


/*! Check if a recursator (i.e. its referred matrix) is empty.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
bool inline is_empty(const recursator<Matrix>& rec)
{
    return rec.is_empty();
}

/*! Check if a recursator (i.e. its referred matrix) fills the
    entire block, i.e. if the number of rows and columns are both
    equal to the virtual bound.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
bool inline is_full(const recursator<Matrix>& rec)
{
    return num_rows(rec) == rec.bound() && num_cols(rec) == rec.bound();
}

/*! The number of rows that a sub-matrix would have if it was constructed.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
typename recursator<Matrix>::size_type
inline num_rows(const recursator<Matrix>& rec)
{
    using std::min;
    typename recursator<Matrix>::size_type tmp= num_rows(rec.my_sub_matrix);
    return rec.my_first_row >= tmp ? 0 : min(rec.my_bound, tmp - rec.my_first_row);
}

/*! The number of columns that a sub-matrix would have if it was constructed.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
typename recursator<Matrix>::size_type
inline num_cols(const recursator<Matrix>& rec)
{
    using std::min;
    typename recursator<Matrix>::size_type tmp= num_cols(rec.my_sub_matrix);
    return rec.my_first_col >= tmp ? 0 : min(rec.my_bound, tmp - rec.my_first_col);
}

/*! The number of elements (rows times columns) that a sub-matrix would have if it was constructed.
    \sa \ref rec_intro "recursion intro"
**/
template <typename Matrix>
typename recursator<Matrix>::size_type
inline size(const recursator<Matrix>& rec)
{
    return num_rows(rec) * num_cols(rec);
}

}} // namespace mtl::matrix

#endif // MTL_MATRIX_RECURATOR_INCLUDE
