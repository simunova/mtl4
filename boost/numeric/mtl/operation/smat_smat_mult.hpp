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

#ifndef MTL_SMAT_SMAT_MULT_INCLUDE
#define MTL_SMAT_SMAT_MULT_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/transposed_matrix_type.hpp>
#include <boost/numeric/mtl/operation/mult_assign_mode.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void smat_smat_mult(const MatrixA& A, const MatrixB& B, MatrixC& C, Assign, 
			   tag::row_major,  // orientation A 
			   tag::row_major)  // orientation B
{
    if (Assign::init_to_zero) set_to_zero(C);
    
    // Average numbers of non-zeros per row
    double ava= num_cols(A) ? double(A.nnz()) / num_cols(A) : 0, 
	   avb= num_rows(B) ? double(B.nnz()) / num_rows(B) : 0; 

    // Define Updater type corresponding to assign mode
    typedef typename Collection<MatrixC>::value_type                            C_value_type;
    typedef typename operations::update_assign_mode<Assign, C_value_type>::type Updater;

    // Reserve 20% over the average's product for entries in C
    mat::inserter<MatrixC, Updater>     ins(C, int( ava * avb * 1.4 ));

    typename traits::row<MatrixA>::type             row_A(A); 
    typename traits::col<MatrixA>::type             col_A(A); 
    typename traits::const_value<MatrixA>::type     value_A(A); 

    typename traits::col<MatrixB>::type             col_B(B); 
    typename traits::const_value<MatrixB>::type     value_B(B); 

    typedef typename traits::range_generator<tag::row, MatrixA>::type  cursor_type;
    cursor_type cursor = begin<tag::row>(A), cend = end<tag::row>(A); 
    for (unsigned ra= 0; cursor != cend; ++ra, ++cursor) {
	// Iterate over non-zeros of each row of A
	typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) {
	    typename Collection<MatrixA>::size_type     ca= col_A(*icursor);   // column of non-zero
	    typename Collection<MatrixA>::value_type    va= value_A(*icursor); // value of non-zero
 
	    // Get cursor corresponding to row 'ca' in matrix B
	    typedef typename traits::range_generator<tag::row, MatrixB>::type  B_cursor_type;
	    B_cursor_type B_cursor = begin<tag::row>(B);
		B_cursor += int(ca); // not elegant but prevents warning

	    // Iterate over non-zeros of this row 
	    typedef typename traits::range_generator<tag::nz, B_cursor_type>::type ib_cursor_type;
	    for (ib_cursor_type ib_cursor = begin<tag::nz>(B_cursor), ib_cend = end<tag::nz>(B_cursor); 
		 ib_cursor != ib_cend; ++ib_cursor) {
		typename Collection<MatrixB>::size_type     cb= col_B(*ib_cursor);   // column of non-zero
		typename Collection<MatrixB>::value_type    vb= value_B(*ib_cursor); // value of non-zero
		ins(ra, cb) << va * vb;		
	    }
	}
    }
}

template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void smat_smat_mult(const MatrixA& A, const MatrixB& B, MatrixC& C, Assign, 
			   tag::col_major,  // orientation A 
			   tag::col_major)  // orientation B
{
    if (Assign::init_to_zero) set_to_zero(C);
    
    // Average numbers of non-zeros per column
    double ava= double(A.nnz()) / num_cols(A), avb= double(B.nnz()) / num_cols(B); 

    // Define Updater type corresponding to assign mode
    typedef typename Collection<MatrixC>::value_type                            C_value_type;
    typedef typename operations::update_assign_mode<Assign, C_value_type>::type Updater;

    // Reserve 20% over the average's product for entries in C
    mat::inserter<MatrixC, Updater>     ins(C, int( ava * avb * 1.2 ));

    typename traits::row<MatrixA>::type             row_A(A); 
    typename traits::col<MatrixA>::type             col_A(A); 
    typename traits::const_value<MatrixA>::type     value_A(A); 

    typename traits::row<MatrixB>::type             row_B(B); 
    typename traits::col<MatrixB>::type             col_B(B); 
    typename traits::const_value<MatrixB>::type     value_B(B); 

    typedef typename traits::range_generator<tag::col, MatrixB>::type  cursor_type;
    cursor_type cursor = begin<tag::col>(B), cend = end<tag::col>(B); 
    for (unsigned cb= 0; cursor != cend; ++cb, ++cursor) {
	// Iterate over non-zeros of each column of B
	typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) {
	    typename Collection<MatrixB>::size_type     rb= row_B(*icursor);   // row of non-zero
	    typename Collection<MatrixB>::value_type    vb= value_B(*icursor); // value of non-zero
 
	    // Get cursor corresponding to column 'rb' in matrix A
	    typedef typename traits::range_generator<tag::col, MatrixA>::type  A_cursor_type;
	    A_cursor_type A_cursor = begin<tag::col>(A);
	    A_cursor+= rb;

	    // Iterate over non-zeros of this column
	    typedef typename traits::range_generator<tag::nz, A_cursor_type>::type ia_cursor_type;
	    for (ia_cursor_type ia_cursor = begin<tag::nz>(A_cursor), ia_cend = end<tag::nz>(A_cursor); 
		 ia_cursor != ia_cend; ++ia_cursor) {
		typename Collection<MatrixA>::size_type     ra= row_A(*ia_cursor);   // row of non-zero
		typename Collection<MatrixA>::value_type    va= value_A(*ia_cursor); // value of non-zero
		ins(ra, cb) << va * vb;		
	    }
	}
    }
}


template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void smat_smat_mult(const MatrixA& A, const MatrixB& B, MatrixC& C, Assign, 
			   tag::col_major,  // orientation A 
			   tag::row_major)  // orientation B
{
    if (Assign::init_to_zero) set_to_zero(C);
    
    // Average numbers of non-zeros per row
    double ava= double(A.nnz()) / num_rows(A), avb= double(B.nnz()) / num_rows(B); 

    // Define Updater type corresponding to assign mode
    typedef typename Collection<MatrixC>::value_type                            C_value_type;
    typedef typename operations::update_assign_mode<Assign, C_value_type>::type Updater;

    // Reserve 20% over the average's product for entries in C
    mat::inserter<MatrixC, Updater>     ins(C, int( ava * avb * 1.2 ));

    typename traits::row<MatrixA>::type             row_A(A); 
    typename traits::col<MatrixA>::type             col_A(A); 
    typename traits::const_value<MatrixA>::type     value_A(A); 

    typename traits::row<MatrixB>::type             row_B(B); 
    typename traits::col<MatrixB>::type             col_B(B); 
    typename traits::const_value<MatrixB>::type     value_B(B); 

    typedef typename traits::range_generator<tag::col, MatrixA>::type  A_cursor_type;
    A_cursor_type A_cursor = begin<tag::col>(A), A_cend = end<tag::col>(A); 

    typedef typename traits::range_generator<tag::row, MatrixB>::type  B_cursor_type;
    B_cursor_type B_cursor = begin<tag::row>(B);

    for (unsigned ca= 0; A_cursor != A_cend; ++ca, ++A_cursor, ++B_cursor) {

	// Iterate over non-zeros of A's column
	typedef typename traits::range_generator<tag::nz, A_cursor_type>::type ia_cursor_type;
	for (ia_cursor_type ia_cursor = begin<tag::nz>(A_cursor), ia_cend = end<tag::nz>(A_cursor); 
	     ia_cursor != ia_cend; ++ia_cursor) 
        {
	    typename Collection<MatrixA>::size_type     ra= row_A(*ia_cursor);   // row of non-zero
	    typename Collection<MatrixA>::value_type    va= value_A(*ia_cursor); // value of non-zero

	    // Iterate over non-zeros of B's row 
	    typedef typename traits::range_generator<tag::nz, B_cursor_type>::type ib_cursor_type;
	    for (ib_cursor_type ib_cursor = begin<tag::nz>(B_cursor), ib_cend = end<tag::nz>(B_cursor); 
		 ib_cursor != ib_cend; ++ib_cursor) 
            {
		typename Collection<MatrixB>::size_type     cb= col_B(*ib_cursor);   // column of non-zero
		typename Collection<MatrixB>::value_type    vb= value_B(*ib_cursor); // value of non-zero
		ins(ra, cb) << va * vb;		
	    }
	}
    }
}


template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign>
inline void smat_smat_mult(const MatrixA& A, const MatrixB& B, MatrixC& C, Assign, 
			   tag::row_major,  // orientation A 
			   tag::col_major)  // orientation B
{
	vampir_trace<4020> tracer;
    // Copy B into a sparse row-major matrix
    typename mtl::traits::transposed_sparse_matrix_type<MatrixB>::type B_copy(num_rows(B), num_cols(B));
    B_copy= B;
    smat_smat_mult(A, B_copy, C, Assign(), tag::row_major(), tag::row_major());
}

} // namespace mtl

#endif // MTL_SMAT_SMAT_MULT_INCLUDE
