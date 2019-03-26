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

#ifndef MTL_SUBMATRIX_INCLUDE
#define MTL_SUBMATRIX_INCLUDE

#include <cmath>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

// Functor type as background for free submatrix function
template <typename Matrix>
struct sub_matrix_t
{
    // typedef  *user_defined*   sub_matrix_type;
    // typedef  *user_defined*   const_sub_matrix_type;
    // typedef  *user_defined*   size_type;
    // sub_matrix_type operator()(Matrix&, size_type, size_type, size_type, size_type);
    // const_sub_matrix_type operator()(Matrix const&, size_type, size_type, size_type, size_type);
};
    
    namespace impl {

	template <typename Matrix>
	inline void correct_sub_matrix_indices(Matrix const& matrix, 
					       typename sub_matrix_t<Matrix>::size_type& begin_row, 
					       typename sub_matrix_t<Matrix>::size_type& end_row, 
					       typename sub_matrix_t<Matrix>::size_type& begin_col, 
					       typename sub_matrix_t<Matrix>::size_type& end_col)
	{
		vampir_trace<3036> tracer;
	    using std::min;
	    MTL_DEBUG_THROW_IF( is_negative(begin_row) || is_negative(end_row), index_out_of_range());
	    end_row= min(end_row, num_rows(matrix));
	    begin_row= min(begin_row, end_row); // implies min(begin_row, num_rows(matrix))
	    
	    MTL_DEBUG_THROW_IF( is_negative(begin_col) || is_negative(end_col), index_out_of_range());
	    end_col= min(end_col, num_cols(matrix));
	    begin_col= min(begin_col, end_col); // implies likewise
	}

    } // namespace impl

///Returns sub-matrix B with begin_row, end_row, begin_col, end_col from %matrix A
template <typename Matrix>
inline typename sub_matrix_t<Matrix>::sub_matrix_type 
sub_matrix(Matrix& matrix, 
	   typename sub_matrix_t<Matrix>::size_type begin_row, 
	   typename sub_matrix_t<Matrix>::size_type end_row, 
	   typename sub_matrix_t<Matrix>::size_type begin_col, 
	   typename sub_matrix_t<Matrix>::size_type end_col)
{
    impl::correct_sub_matrix_indices(matrix, begin_row, end_row, begin_col, end_col);
    return sub_matrix_t<Matrix>()(matrix, begin_row, end_row, begin_col, end_col);
}

template <typename Matrix>
inline typename sub_matrix_t<Matrix>::const_sub_matrix_type 
sub_matrix(Matrix const& matrix, 
	   typename sub_matrix_t<Matrix>::size_type begin_row, 
	   typename sub_matrix_t<Matrix>::size_type end_row, 
	   typename sub_matrix_t<Matrix>::size_type begin_col, 
	   typename sub_matrix_t<Matrix>::size_type end_col)
{
    impl::correct_sub_matrix_indices(matrix, begin_row, end_row, begin_col, end_col);
    return sub_matrix_t<Matrix>()(matrix, begin_row, end_row, begin_col, end_col);
}

}} // namespace mtl::matrix

#endif // MTL_SUBMATRIX_INCLUDE
