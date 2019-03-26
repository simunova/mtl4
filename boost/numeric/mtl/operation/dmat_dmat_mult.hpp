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

#ifndef MTL_DMAT_DMAT_MULT_INCLUDE
#define MTL_DMAT_DMAT_MULT_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/operation/cursor_pseudo_dot.hpp>
#include <boost/numeric/mtl/operation/multi_action_block.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/operation/static_size.hpp>
#include <boost/numeric/mtl/operation/static_num_rows.hpp>
#include <boost/numeric/mtl/operation/static_num_cols.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/flatcat.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/meta_math/loop.hpp>
#include <boost/numeric/mtl/recursion/base_case_test.hpp>
#include <boost/numeric/mtl/recursion/base_case_matrix.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>
#include <boost/numeric/mtl/recursion/base_case_cast.hpp>
#include <boost/numeric/mtl/interface/blas.hpp>

#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/no_op.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#include <iostream>
#include <complex>

namespace mtl {

// =====================================================
// Generic matrix product with cursors and property maps
// =====================================================
   

// To allow 5th parameter, is ignored
// This is the bottom line of dmat_dmat_mult implementations.
// All MTL4 matrix types and views (so far) have cursors and property maps
// so that we disabled the backup functor by default.
// If some type has no cursor, one can still use a backup functor (whatever this may be).
template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign= assign::assign_sum,
	  typename Backup= no_op> 
struct gen_cursor_dmat_dmat_mult_ft
{
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	apply(A, B, C, traits::flatcat1<MatrixA, tag::has_cursor>(), traits::flatcat1<MatrixB, tag::has_cursor>());
    }   

private:
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::universe, tag::universe)
    {
	Backup()(A, B, C);
    }

    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::flat<tag::has_cursor>, tag::flat<tag::has_cursor>)
    {
	// asm("#cursor");
	vampir_trace<4001> tracer;
	// std::cout << "Canonical cursor\n";
	typedef glas::tag::row                                          row;
	typedef glas::tag::col                                          col;
	typedef glas::tag::all                                          all;

        typedef typename traits::const_value<MatrixA>::type                a_value_type;
        typedef typename traits::const_value<MatrixB>::type                b_value_type;
        typedef typename traits::value<MatrixC>::type                      c_value_type;

        typedef typename traits::range_generator<row, MatrixA>::type     a_cur_type;
        typedef typename traits::range_generator<row, MatrixC>::type     c_cur_type;
        
        typedef typename traits::range_generator<col, MatrixB>::type     b_cur_type;
        typedef typename traits::range_generator<all, c_cur_type>::type  c_icur_type;

        typedef typename traits::range_generator<all, a_cur_type>::type  a_icur_type;
        typedef typename traits::range_generator<all, b_cur_type>::type  b_icur_type;

#ifndef NDEBUG
	typename traits::row<MatrixA>::type                              row_a(A); 
	typename traits::col<MatrixA>::type                              col_a(A); 
	typename traits::row<MatrixB>::type                              row_b(B); 
	typename traits::col<MatrixB>::type                              col_b(B); 
	typename traits::row<MatrixC>::type                              row_c(C); 
	typename traits::col<MatrixC>::type                              col_c(C); 
#else
#  undef MTL_DEBUG_DMAT_DMAT_MULT // doesn't work with NDEBUG
#endif

	if (Assign::init_to_zero) set_to_zero(C);

	a_value_type   a_value(A);
	b_value_type   b_value(B);
	c_value_type   c_value(C);
    		
	a_cur_type ac= begin<row>(A), aend= end<row>(A);
	for (c_cur_type cc= begin<row>(C); ac != aend; ++ac, ++cc) {
	    
	    b_cur_type bc= begin<col>(B), bend= end<col>(B);
	    for (c_icur_type cic= begin<all>(cc); bc != bend; ++bc, ++cic) { 
		
		typename MatrixC::value_type c_tmp(c_value(*cic));
#ifdef MTL_DEBUG_DMAT_DMAT_MULT
		std::cout << "Calculating C[" << row_c(*cic) << "][" << col_c(*cic) << "], initial value is "
			  << c_tmp << "\n";
#endif
		a_icur_type aic= begin<all>(ac), aiend= end<all>(ac); 
		for (b_icur_type bic= begin<all>(bc); aic != aiend; ++aic, ++bic) {
#ifdef MTL_DEBUG_DMAT_DMAT_MULT
		    std::cout << "Updating with A[" << row_a(*aic) << "][" << col_a(*aic) << "] /* value is "
			      << a_value(*aic) << " */ * ";
		    std::cout << "B[" << row_b(*bic) << "][" << col_b(*bic) << "] /* value is "
			      << b_value(*bic) << " */ * ";
#endif
		    assert(row_a(*aic) == row_c(*cic)); // Must do here because ac has no props
		    assert(col_b(*bic) == col_c(*cic));
		    assert(col_a(*aic) == row_b(*bic));
		    Assign::update(c_tmp, a_value(*aic) * b_value(*bic));
#ifdef MTL_DEBUG_DMAT_DMAT_MULT
		    std::cout << "C's current value is " << c_tmp << "\n";
#endif
		}
		c_value(*cic, c_tmp);
	    }
	} 
    }
};


template <typename Assign= assign::assign_sum,
	  typename Backup= no_op>     // To allow 2nd parameter, is ignored
struct gen_cursor_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_cursor_dmat_dmat_mult_ft<MatrixA, MatrixB, MatrixC, Assign, Backup>()(A, B, C);
    }
};


// =====================================
// Generic matrix product with iterators
// =====================================

template <typename MatrixA, typename MatrixB, typename MatrixC, 
	  typename Assign= assign::assign_sum, 
	  typename Backup= gen_cursor_dmat_dmat_mult_t<Assign> > 
struct gen_dmat_dmat_mult_ft
{
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	apply(A, B, C, traits::flatcat1<MatrixA, tag::has_iterator>(), traits::flatcat1<MatrixB, tag::has_iterator>());
    }   

private:
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::universe, tag::universe)
    {
	Backup()(A, B, C);
    }

    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::flat<tag::has_iterator>, tag::flat<tag::has_iterator>)
    {
	// asm("#iterator");
	vampir_trace<4002> tracer;
	// std::cout << "Canonical iterator\n";
	using namespace tag;
	using traits::range_generator;  
        typedef typename range_generator<row, MatrixA>::type       a_cur_type;             
        typedef typename range_generator<row, MatrixC>::type       c_cur_type;             
	typedef typename range_generator<col, MatrixB>::type       b_cur_type;             
        typedef typename range_generator<iter::all, c_cur_type>::type   c_icur_type;            
        typedef typename range_generator<const_iter::all, a_cur_type>::type  a_icur_type;            
        typedef typename range_generator<const_iter::all, b_cur_type>::type  b_icur_type;          

	if (Assign::init_to_zero) set_to_zero(C);

	a_cur_type ac= mtl::begin<row>(A), aend= mtl::end<row>(A);
	for (c_cur_type cc= mtl::begin<row>(C); ac != aend; ++ac, ++cc) {

	    b_cur_type bc= mtl::begin<col>(B), bend= mtl::end<col>(B);
	    for (c_icur_type cic= mtl::begin<iter::all>(cc); bc != bend; ++bc, ++cic) { 
		    
		typename MatrixC::value_type c_tmp(*cic);
		a_icur_type aic= mtl::begin<const_iter::all>(ac), aiend= mtl::end<const_iter::all>(ac); 
		for (b_icur_type bic= mtl::begin<const_iter::all>(bc); aic != aiend; ++aic, ++bic) {
		    Assign::update(c_tmp, *aic * *bic);
		}
		*cic= c_tmp;
	    }
	}
    }    
};


template <typename Assign= assign::assign_sum,
	  typename Backup= gen_cursor_dmat_dmat_mult_t<Assign> >   
struct gen_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_dmat_dmat_mult_ft<MatrixA, MatrixB, MatrixC, Assign, Backup>()(A, B, C);
    }
};


/*

Unrolling matrix product with dimensions that are not multiples of blocks

1. Do with optimization:
   C_nw += A_nw * B_nw
   - wherby the matrix dimensions of sub-matrices are the largest multiples of block sizes 
     smaller or equal to the matrix dimensions of the original matrix


2. Do without optimization
   C_nw += A_ne * B_sw
   C_ne += A_n * B_e
   C_s += A_s * B

The inner loop can be unrolled arbitrarily. So, we can simplify

1. Do with optimization:
   C_nw += A_n * B_w
   - wherby the matrix dimensions of sub-matrices are the largest multiples of block sizes 
     smaller or equal to the matrix dimensions of the original matrix


2. Do with optimization only in inner loop
   C_ne += A_n * B_e
   C_s += A_s * B
  

*/

// =======================
// Unrolled with iterators
// required has_2D_layout
// =======================

// Define defaults if not yet given as Compiler flag
#ifndef MTL_DMAT_DMAT_MULT_TILING1
#  define MTL_DMAT_DMAT_MULT_TILING1 2
#endif

#ifndef MTL_DMAT_DMAT_MULT_TILING2
#  define MTL_DMAT_DMAT_MULT_TILING2 4
#endif

#ifndef MTL_DMAT_DMAT_MULT_INNER_UNROLL
#  define MTL_DMAT_DMAT_MULT_INNER_UNROLL 8
#endif


template <unsigned long Index0, unsigned long Max0, unsigned long Index1, unsigned long Max1, typename Assign>
struct gen_tiling_dmat_dmat_mult_block
  : public meta_math::loop2<Index0, Max0, Index1, Max1>
{
    typedef meta_math::loop2<Index0, Max0, Index1, Max1>                              base;
    typedef gen_tiling_dmat_dmat_mult_block<base::next_index0, Max0, base::next_index1, Max1, Assign>  next_t;

    template <typename Value, typename ValueA, typename SizeA, typename ValueB, typename SizeB>
    static inline void apply(Value& tmp00, Value& tmp01, Value& tmp02, Value& tmp03, Value& tmp04, 
			     Value& tmp05, Value& tmp06, Value& tmp07, Value& tmp08, Value& tmp09, 
			     Value& tmp10, Value& tmp11, Value& tmp12, Value& tmp13, Value& tmp14, Value& tmp15, 
			     ValueA *begin_a, SizeA& ari, ValueB *begin_b, SizeB& bci)
    {
	tmp00+= begin_a[ base::index0 * ari ] * begin_b[ base::index1 * bci ];
	next_t::apply(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		      tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
		      begin_a, ari, begin_b, bci); 
    }

    template <typename Value, typename MatrixC, typename SizeC>
    static inline void update(Value& tmp00, Value& tmp01, Value& tmp02, Value& tmp03, Value& tmp04, 
			      Value& tmp05, Value& tmp06, Value& tmp07, Value& tmp08, Value& tmp09, 
			      Value& tmp10, Value& tmp11, Value& tmp12, Value& tmp13, Value& tmp14, Value& tmp15,
			      MatrixC& C, SizeC i, SizeC k)
    {
	Assign::update(C(i + base::index0, k + base::index1), tmp00);
	next_t::update(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		       tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
		       C, i, k);
    }
};

template <unsigned long Max0, unsigned long Max1, typename Assign>
struct gen_tiling_dmat_dmat_mult_block<Max0, Max0, Max1, Max1, Assign>
    : public meta_math::loop2<Max0, Max0, Max1, Max1>
{
    typedef meta_math::loop2<Max0, Max0, Max1, Max1>  base;

    template <typename Value, typename ValueA, typename SizeA, typename ValueB, typename SizeB>
    static inline void apply(Value& tmp00, Value&, Value&, Value&, Value&, 
			     Value&, Value&, Value&, Value&, Value&, 
			     Value&, Value&, Value&, Value&, Value&, Value&, 
			     ValueA *begin_a, SizeA& ari, ValueB *begin_b, SizeB& bci)
    {
	tmp00+= begin_a[ base::index0 * ari ] * begin_b[ base::index1 * bci ];
    }

    template <typename Value, typename MatrixC, typename SizeC>
    static inline void update(Value& tmp00, Value&, Value&, Value&, Value&, 
			      Value&, Value&, Value&, Value&, Value&, 
			      Value&, Value&, Value&, Value&, Value&, Value&,
			      MatrixC& C, SizeC i, SizeC k)
    {
	Assign::update(C(i + base::index0, k + base::index1), tmp00);
    }
};


template <typename MatrixA, typename MatrixB, typename MatrixC,
	  unsigned long Tiling1= MTL_DMAT_DMAT_MULT_TILING1,
	  unsigned long Tiling2= MTL_DMAT_DMAT_MULT_TILING2,
	  typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_tiling_dmat_dmat_mult_ft
{
    MTL_STATIC_ASSERT((Tiling1 * Tiling2 <= 16), "Tile (Tiling1 * Tiling2) cannot be larger than 16.");
  
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	apply(A, B, C, traits::layout_flatcat<MatrixA>(), traits::layout_flatcat<MatrixB>());
    }   
 
private:
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::universe, tag::universe)
    {
	Backup()(A, B, C);
    }

    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::flat<tag::has_2D_layout>, tag::flat<tag::has_2D_layout>)
    {
	// asm("#tiling");
	vampir_trace<4003> tracer;
	// Indices run out of range for smaller matrices
	if (num_rows(A) < 2 || num_cols(A) < 2 || num_cols(B) < 2) {
	    Backup()(A, B, C);
	    return;
	}

	// std::cout << "meta-unrolling\n";
	if (Assign::init_to_zero) set_to_zero(C);

	typedef gen_tiling_dmat_dmat_mult_block<1, Tiling1, 1, Tiling2, Assign>  block;
	typedef typename MatrixC::size_type                                          size_type;
	typedef typename MatrixC::value_type                                         value_type;
	const value_type z= math::zero(C[0][0]);    // if this are matrices we need their size

	size_type i_max= num_rows(C), i_block= Tiling1 * (i_max / Tiling1),
	          k_max= num_cols(C), k_block= Tiling2 * (k_max / Tiling2);
	size_t ari= &A(1, 0) - &A(0, 0), // how much is the offset of A's entry increased by incrementing row
	       aci= &A(0, 1) - &A(0, 0), bri= &B(1, 0) - &B(0, 0), bci= &B(0, 1) - &B(0, 0);
	    
	// C_nw += A_nw * B_nw
	for (size_type i= 0; i < i_block; i+= Tiling1)
	    for (size_type k= 0; k < k_block; k+= Tiling2) {

		value_type tmp00= z, tmp01= z, tmp02= z, tmp03= z, tmp04= z,
                           tmp05= z, tmp06= z, tmp07= z, tmp08= z, tmp09= z,
 		           tmp10= z, tmp11= z, tmp12= z, tmp13= z, tmp14= z, tmp15= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    block::apply(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
				 tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
				 begin_a, ari, begin_b, bci); 
		block::update(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
			      tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
			      C, i, k);
	    }

	// C_ne += A_n * B_e
	for (size_type i= 0; i < i_block; i++)
	    for (size_type k = k_block; k < k_max; k++) {
		value_type tmp00= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    tmp00 += *begin_a * *begin_b;
		Assign::update(C(i, k), tmp00);
	    }

	// C_s += A_s * B
	for (size_type i= i_block; i < i_max; i++)
	    for (size_type k = 0; k < k_max; k++) {
		value_type tmp00= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    tmp00 += *begin_a * *begin_b;
		Assign::update(C(i, k), tmp00);
	    }
    }
};

template <unsigned long Tiling1= MTL_DMAT_DMAT_MULT_TILING1,
	  unsigned long Tiling2= MTL_DMAT_DMAT_MULT_TILING2,
	  typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_tiling_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_tiling_dmat_dmat_mult_ft<
	     MatrixA, MatrixB, MatrixC, Tiling1, Tiling2, Assign, Backup
	>()(A, B, C);
    }
};



// =================================
// Unrolled with iterators fixed 4x4
// required has_2D_layout
// =================================


template <typename MatrixA, typename MatrixB, typename MatrixC, 
	  typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_tiling_44_dmat_dmat_mult_ft
{
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	apply(A, B, C, traits::layout_flatcat<MatrixA>(), traits::layout_flatcat<MatrixB>());
    }   
 
private:
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::universe, tag::universe)
    {
	Backup()(A, B, C);
    }

    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::flat<tag::has_2D_layout>, tag::flat<tag::has_2D_layout>)
    {
	// asm("#tiling44");
	vampir_trace<4004> tracer;
	// Indices run out of range for smaller matrices
	if (num_rows(A) < 2 || num_cols(A) < 2 || num_cols(B) < 2) {
	    Backup()(A, B, C);
	    return;
	}

        // std::cout << "4x4 unrolling\n";
	if (Assign::init_to_zero) set_to_zero(C);

	typedef typename MatrixC::size_type                                          size_type;
	typedef typename MatrixC::value_type                                         value_type;

	const size_type  Tiling1= 4, Tiling2= 4;
	const value_type z= math::zero(C[0][0]);    // if this are matrices we need their size

	size_type i_max= num_rows(C), i_block= Tiling1 * (i_max / Tiling1),
	          k_max= num_cols(C), k_block= Tiling2 * (k_max / Tiling2);
	size_t ari= &A(1, 0) - &A(0, 0), // how much is the offset of A's entry increased by incrementing row
	       aci= &A(0, 1) - &A(0, 0), bri= &B(1, 0) - &B(0, 0), bci= &B(0, 1) - &B(0, 0);

	// C_nw += A_nw * B_nw
	for (size_type i= 0; i < i_block; i+= Tiling1)
	    for (size_type k= 0; k < k_block; k+= Tiling2) {

		value_type tmp00= z, tmp01= z, tmp02= z, tmp03= z, tmp04= z,
                           tmp05= z, tmp06= z, tmp07= z, tmp08= z, tmp09= z,
 		           tmp10= z, tmp11= z, tmp12= z, tmp13= z, tmp14= z, tmp15= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri) {
		    tmp00+= begin_a[ 0 * ari ] * begin_b[ 0 * bci ];
		    tmp01+= begin_a[ 0 * ari ] * begin_b[ 1 * bci ];
		    tmp02+= begin_a[ 0 * ari ] * begin_b[ 2 * bci ];
		    tmp03+= begin_a[ 0 * ari ] * begin_b[ 3 * bci ];
		    tmp04+= begin_a[ 1 * ari ] * begin_b[ 0 * bci ];
		    tmp05+= begin_a[ 1 * ari ] * begin_b[ 1 * bci ];
		    tmp06+= begin_a[ 1 * ari ] * begin_b[ 2 * bci ];
		    tmp07+= begin_a[ 1 * ari ] * begin_b[ 3 * bci ];
		    tmp08+= begin_a[ 2 * ari ] * begin_b[ 0 * bci ];
		    tmp09+= begin_a[ 2 * ari ] * begin_b[ 1 * bci ];
		    tmp10+= begin_a[ 2 * ari ] * begin_b[ 2 * bci ];
		    tmp11+= begin_a[ 2 * ari ] * begin_b[ 3 * bci ];
		    tmp12+= begin_a[ 3 * ari ] * begin_b[ 0 * bci ];
		    tmp13+= begin_a[ 3 * ari ] * begin_b[ 1 * bci ];
		    tmp14+= begin_a[ 3 * ari ] * begin_b[ 2 * bci ];
		    tmp15+= begin_a[ 3 * ari ] * begin_b[ 3 * bci ];
		}
		Assign::update(C(i + 0, k + 0), tmp00);
		Assign::update(C(i + 0, k + 1), tmp01);
		Assign::update(C(i + 0, k + 2), tmp02);
		Assign::update(C(i + 0, k + 3), tmp03);
		Assign::update(C(i + 1, k + 0), tmp04);
		Assign::update(C(i + 1, k + 1), tmp05);
		Assign::update(C(i + 1, k + 2), tmp06);
		Assign::update(C(i + 1, k + 3), tmp07);
		Assign::update(C(i + 2, k + 0), tmp08);
		Assign::update(C(i + 2, k + 1), tmp09);
		Assign::update(C(i + 2, k + 2), tmp10);
		Assign::update(C(i + 2, k + 3), tmp11);
		Assign::update(C(i + 3, k + 0), tmp12);
		Assign::update(C(i + 3, k + 1), tmp13);
		Assign::update(C(i + 3, k + 2), tmp14);
		Assign::update(C(i + 3, k + 3), tmp15);
	    }

	// C_ne += A_n * B_e
	for (size_type i= 0; i < i_block; i++)
	    for (size_type k = k_block; k < k_max; k++) {
		value_type tmp00= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    tmp00 += *begin_a * *begin_b;
		Assign::update(C(i, k), tmp00);
	    }

	// C_s += A_s * B
	for (size_type i= i_block; i < i_max; i++)
	    for (size_type k = 0; k < k_max; k++) {
		value_type tmp00= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    tmp00 += *begin_a * *begin_b;
		Assign::update(C(i, k), tmp00);
	    }
    }
};

template <typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_tiling_44_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_tiling_44_dmat_dmat_mult_ft<
	     MatrixA, MatrixB, MatrixC, Assign, Backup
	>()(A, B, C);
    }
};




// =================================
// Unrolled with iterators fixed 2x2
// required has_2D_layout
// =================================


template <typename MatrixA, typename MatrixB, typename MatrixC, 
	  typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_tiling_22_dmat_dmat_mult_ft
{
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	apply(A, B, C,  traits::layout_flatcat<MatrixA>(), traits::layout_flatcat<MatrixB>());
    }   
 
private:
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::universe, tag::universe)
    {
	Backup()(A, B, C);
    }

    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::flat<tag::has_2D_layout>, tag::flat<tag::has_2D_layout>)
    {
	// asm("#tiling22");
	vampir_trace<4005> tracer;
	// Indices run out of range for smaller matrices
	if (num_rows(A) < 2 || num_cols(A) < 2 || num_cols(B) < 2) {
	    Backup()(A, B, C);
	    return;
	}

        // std::cout << "2x2 unrolling\n";
	if (Assign::init_to_zero) set_to_zero(C);

	typedef typename MatrixC::size_type                                          size_type;
	typedef typename MatrixC::value_type                                         value_type;

	const size_type  Tiling1= 2, Tiling2= 2;
	const value_type z= math::zero(C[0][0]);    // if this are matrices we need their size

	size_type i_max= num_rows(C), i_block= Tiling1 * (i_max / Tiling1),
	          k_max= num_cols(C), k_block= Tiling2 * (k_max / Tiling2);
	size_t ari= &A(1, 0) - &A(0, 0), // how much is the offset of A's entry increased by incrementing row
	       aci= &A(0, 1) - &A(0, 0), bri= &B(1, 0) - &B(0, 0), bci= &B(0, 1) - &B(0, 0);

	// C_nw += A_nw * B_nw
	for (size_type i= 0; i < i_block; i+= Tiling1)
	    for (size_type k= 0; k < k_block; k+= Tiling2) {

		value_type tmp00= z, tmp01= z, tmp02= z, tmp03= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri) {
		    tmp00+= begin_a[ 0 ] * begin_b[ 0 ];
		    tmp01+= begin_a[ 0 ] * begin_b[bci];
		    tmp02+= begin_a[ari] * begin_b[ 0 ];
		    tmp03+= begin_a[ari] * begin_b[bci];
		}
		Assign::update(C(i + 0, k + 0), tmp00);
		Assign::update(C(i + 0, k + 1), tmp01);
		Assign::update(C(i + 1, k + 0), tmp02);
		Assign::update(C(i + 1, k + 1), tmp03);
	    }

	// C_ne += A_n * B_e
	for (size_type i= 0; i < i_block; i++)
	    for (size_type k = k_block; k < k_max; k++) {
		value_type tmp00= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    tmp00 += *begin_a * *begin_b;
		Assign::update(C(i, k), tmp00);
	    }

	// C_s += A_s * B
	for (size_type i= i_block; i < i_max; i++)
	    for (size_type k = 0; k < k_max; k++) {
		value_type tmp00= z;
		const typename MatrixA::value_type *begin_a= &A(i, 0), *end_a= begin_a + num_cols(A) * aci;
		const typename MatrixB::value_type *begin_b= &B(0, k);

		for (; begin_a != end_a; begin_a+= aci, begin_b+= bri)
		    tmp00 += *begin_a * *begin_b;
		Assign::update(C(i, k), tmp00);
	    }
    }
};

template <typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_tiling_22_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_tiling_22_dmat_dmat_mult_ft<
	     MatrixA, MatrixB, MatrixC, Assign, Backup
	>()(A, B, C);
    }
};




// ========================
// Recursive Multiplication
// ========================

namespace wrec {

    template <typename BaseMult, typename BaseTest= recursion::bound_test_static<64> >
    struct gen_dmat_dmat_mult_t
    {
	template <typename RecA, typename RecB, typename RecC>
	void operator()(RecA const& rec_a, RecB const& rec_b, RecC& rec_c)
	{
	    vampir_trace<4006> tracer;
	    // std::cout << "wrec::mult \n";
	    using namespace recursion;
	    // using mtl::mat::is_empty; // ambiguity with std::tr1::is_empty in VS2010
	    // Ambiguity with boost::is_empty in Open64
	    if (mtl::mat::is_empty(rec_a) || mtl::mat::is_empty(rec_b) || mtl::mat::is_empty(rec_c))
		return;

	    if (BaseTest()(rec_a)) {
		//std::cout << "base_case: A =\n" << *rec_a << "B =\n" << *rec_b;
#               if 0 // ndef NDEBUG
		    typename RecC::sub_matrix_type C_tmp(*rec_c);
		    typename base_case_matrix<typename RecC::sub_matrix_type, BaseTest>::type
			C = base_case_cast<BaseTest>(C_tmp);
#               else
    		    typename base_case_matrix<typename RecC::matrix_type, BaseTest>::type
			C= base_case_cast<BaseTest>(*rec_c);
#               endif
		BaseMult()(base_case_cast<BaseTest>(*rec_a),
			   base_case_cast<BaseTest>(*rec_b), C);
		//std::cout << "base_case: C =\n" << C << "*rec_c = \n" << *rec_c;
#               if 0 //ndef NDEBUG
		    // MTL_ASSERT((C(0, 0) == (*rec_c)(0, 0)), "Result was computed in a copy not in a view/reference!"); // does not compile in set_to_zero_with_user_value_test
#               endif
	    } else {
		RecC c_north_west= north_west(rec_c), c_north_east= north_east(rec_c),
		     c_south_west= south_west(rec_c), c_south_east= south_east(rec_c);

		(*this)(north_west(rec_a), north_west(rec_b), c_north_west);
		(*this)(north_west(rec_a), north_east(rec_b), c_north_east);
		(*this)(south_west(rec_a), north_east(rec_b), c_south_east);
		(*this)(south_west(rec_a), north_west(rec_b), c_south_west);
		(*this)(south_east(rec_a), south_west(rec_b), c_south_west);
		(*this)(south_east(rec_a), south_east(rec_b), c_south_east);
		(*this)(north_east(rec_a), south_east(rec_b), c_north_east);
		(*this)(north_east(rec_a), south_west(rec_b), c_north_west);
	    }
	}
    };

} // namespace wrec


template <typename BaseMult, 
	  typename BaseTest= recursion::bound_test_static<64>,
	  typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_recursive_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	apply(A, B, C, traits::flatcat1<MatrixA, tag::qsub_divisible>(), 
	      traits::flatcat1<MatrixB, tag::qsub_divisible>(), traits::flatcat1<MatrixC, tag::qsub_divisible>());
    }   
 
private:
    // If one matrix is not sub-divisible then take backup function
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, tag::universe, tag::universe, tag::universe)
    {
	Backup()(A, B, C);
    }

    // Only if matrix is sub-divisible, otherwise backup
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, 
	       tag::flat<tag::qsub_divisible>, tag::flat<tag::qsub_divisible>, tag::flat<tag::qsub_divisible>)
    {
	vampir_trace<4007> tracer;
	// std::cout << "do recursion\n";
	if (Assign::init_to_zero) set_to_zero(C);

	// Make sure that mult functor of basecase has appropriate assign mode (in all nestings)
	// i.e. replace assign::assign_sum by assign::plus_sum including backup functor
	
	using mat::recursator;
	recursator<MatrixA>    rec_a(A);
	recursator<MatrixB>    rec_b(B);
	recursator<MatrixC>    rec_c(C);
	equalize_depth(rec_a, rec_b, rec_c);
	
	wrec::gen_dmat_dmat_mult_t<BaseMult, BaseTest>() (rec_a, rec_b, rec_c);
    }
};



// ==================================
// Plattform specific implementations
// ==================================

// Here only general definition that calls backup function
// Special implementations needed in other files, which are included at the end

template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_platform_dmat_dmat_mult_ft
    : public Backup
{};


template <typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_platform_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_platform_dmat_dmat_mult_ft<MatrixA, MatrixB, MatrixC, Assign, Backup>()(A, B, C);
    }
};


// ==================================
// BLAS functions as far as supported
// ==================================


template <typename MatrixA, typename MatrixB, typename MatrixC, typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_blas_dmat_dmat_mult_ft
    : public Backup
{};


#ifdef MTL_HAS_BLAS

namespace detail {

    // Transform from assign representation to BLAS
    double inline dgemm_alpha(assign::assign_sum) { return 1.0; }
    double inline dgemm_alpha(assign::plus_sum) { return 1.0; } 
    double inline dgemm_alpha(assign::minus_sum) { return -1.0; }

    // Transform from assign representation to BLAS
    double inline dgemm_beta(assign::assign_sum) { return 0.0; }
    double inline dgemm_beta(assign::plus_sum) { return 1.0; }
    double inline dgemm_beta(assign::minus_sum) { return 1.0; }

    template <typename Value, typename ParaA, typename ParaB, typename ParaC, typename Function, typename Assign>
    void inline xgemm(const dense2D<Value, ParaA>& A, const dense2D<Value, ParaB>& B, 
		      dense2D<Value, ParaC>& C, Function f, Assign)
    {
	vampir_trace<4008> tracer;
	// std::cout << "use generic BLAS\n";
	int m= num_rows(A), n= num_cols(B), k= num_cols(A), lda= A.get_ldim(), ldb= B.get_ldim(), ldc= C.get_ldim();
	Value alpha= dgemm_alpha(Assign()), beta= dgemm_beta(Assign());
	char a_trans= traits::is_row_major<ParaA>::value ? 'T' : 'N', 
             b_trans= traits::is_row_major<ParaB>::value ? 'T' : 'N';

	if (traits::is_row_major<ParaC>::value) {
	    // C^T= B^T * A^T
	    a_trans= 'T' + 'N' - a_trans; b_trans= 'T' + 'N' - b_trans; 
	    f(&b_trans, &a_trans, &n /* col(B) */, &m /* row(A) */, &k /* col(A)=row(B) */, 
				 &alpha, &B[0][0], &ldb, &A[0][0], &lda, &beta, &C[0][0], &ldc);
	} else 
	    f(&a_trans, &b_trans, &m, &n, &k, &alpha, &A[0][0], &lda, &B[0][0], &ldb, &beta, &C[0][0], &ldc);
    }

} // detail
 
template<typename ParaA, typename ParaB, typename ParaC, typename Assign, typename Backup>
struct gen_blas_dmat_dmat_mult_ft<dense2D<float, ParaA>, dense2D<float, ParaB>, 
				  dense2D<float, ParaC>, Assign, Backup>
{
    void operator()(const dense2D<float, ParaA>& A, const dense2D<float, ParaB>& B, 
		    dense2D<float, ParaC>& C)
    {
	detail::xgemm(A, B, C, MTL_BLAS_NAME(sgemm), Assign());
    }
};

template<typename ParaA, typename ParaB, typename ParaC, typename Assign, typename Backup>
struct gen_blas_dmat_dmat_mult_ft<dense2D<double, ParaA>, dense2D<double, ParaB>, 
				     dense2D<double, ParaC>, Assign, Backup>
{
    void operator()(const dense2D<double, ParaA>& A, const dense2D<double, ParaB>& B, 
		    dense2D<double, ParaC>& C)
    {
	detail::xgemm(A, B, C, MTL_BLAS_NAME(dgemm), Assign());
    }
};


template<typename ParaA, typename ParaB, typename ParaC, typename Assign, typename Backup>
struct gen_blas_dmat_dmat_mult_ft<dense2D<std::complex<float>, ParaA>, dense2D<std::complex<float>, ParaB>, 
				  dense2D<std::complex<float>, ParaC>, Assign, Backup>
{
    void operator()(const dense2D<std::complex<float>, ParaA>& A, const dense2D<std::complex<float>, ParaB>& B, 
		    dense2D<std::complex<float>, ParaC>& C)
    {
	detail::xgemm(A, B, C, MTL_BLAS_NAME(cgemm), Assign());
    }
};


template<typename ParaA, typename ParaB, typename ParaC, typename Assign, typename Backup>
struct gen_blas_dmat_dmat_mult_ft<dense2D<std::complex<double>, ParaA>, dense2D<std::complex<double>, ParaB>, 
				  dense2D<std::complex<double>, ParaC>, Assign, Backup>
{
    void operator()(const dense2D<std::complex<double>, ParaA>& A, const dense2D<std::complex<double>, ParaB>& B, 
		    dense2D<std::complex<double>, ParaC>& C)
    {
	detail::xgemm(A, B, C, MTL_BLAS_NAME(zgemm), Assign());
    }
};



#endif // MTL_HAS_BLAS

template <typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct gen_blas_dmat_dmat_mult_t
  : public Backup
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	gen_blas_dmat_dmat_mult_ft<MatrixA, MatrixB, MatrixC, Assign, Backup>()(A, B, C);
    }
};


// ============================================================
// Switch between two functors depending on minimal matrix size
// ============================================================

template <std::size_t SizeLimit, typename FunctorSmall, typename FunctorLarge>
struct size_switch_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	const bool all_static= traits::is_static<MatrixA>::value && traits::is_static<MatrixB>::value 
	                       && traits::is_static<MatrixC>::value;
	apply(A, B, C, boost::mpl::bool_<all_static>());
    }
    
  private:
    // Decided at run time
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, boost::mpl::false_)
    {
	if (mtl::mat::size(A) <= SizeLimit || mtl::mat::size(B) <= SizeLimit || mtl::mat::size(C) <= SizeLimit)
	    FunctorSmall()(A, B, C);
	else
	    FunctorLarge()(A, B, C);
    }

    // Decided at compile time
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& A, MatrixB const& B, MatrixC& C, boost::mpl::true_)
    {
	const bool is_small= mtl::static_size<MatrixA>::value <= SizeLimit || mtl::static_size<MatrixB>::value <= SizeLimit 
	    || mtl::static_size<MatrixC>::value <= SizeLimit;
	typename boost::mpl::if_c<is_small, FunctorSmall, FunctorLarge>::type()(A, B, C);
    }
};


// ====================================================================================
// Switch between two functors depending on whether matrix size is know at compile time
// ====================================================================================

template <bool IsStatic, typename FunctorStatic, typename FunctorDynamic>
struct static_switch_dmat_dmat_mult_t
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	typename boost::mpl::if_c<IsStatic, FunctorStatic, FunctorDynamic>::type()(A, B, C);
    }
};

// ==============================================================
// Completely unroll fixed size computations, will be tuned later
// ==============================================================

template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1, typename Assign>
struct fully_unroll_dmat_dmat_mult_init_block
  : public meta_math::loop2<Index0, Max0, Index1, Max1>
{
    typedef meta_math::loop2<Index0, Max0, Index1, Max1>                              base;
    typedef fully_unroll_dmat_dmat_mult_init_block<base::next_index0, Max0, base::next_index1, Max1, Assign>  next_t;

    template <typename MatrixA, typename MatrixB, typename MatrixC>
    static inline void apply(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	// Assign::first_update(C[base::index0][base::index1], A[base::index0][0] * B[0][base::index1]);
	Assign::first_update(C(base::index0, base::index1), A(base::index0, 0) * B(0, base::index1));
	next_t::apply(A, B, C);
    }
};

template <std::size_t Max0, std::size_t Max1, typename Assign>
struct fully_unroll_dmat_dmat_mult_init_block<Max0, Max0, Max1, Max1, Assign>
  : public meta_math::loop2<Max0, Max0, Max1, Max1>
{
    typedef meta_math::loop2<Max0, Max0, Max1, Max1>                              base;

    template <typename MatrixA, typename MatrixB, typename MatrixC>
    static inline void apply(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	Assign::first_update(C(base::index0, base::index1), A(base::index0, 0) * B(0, base::index1));
    }
};


template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1, 
	  std::size_t Index2, std::size_t Max2, typename Assign>
struct fully_unroll_dmat_dmat_mult_block
  : public meta_math::loop3<Index0, Max0, Index1, Max1, Index2, Max2>
{
    typedef meta_math::loop3<Index0, Max0, Index1, Max1, Index2, Max2>                              base;
    typedef fully_unroll_dmat_dmat_mult_block<base::next_index0, Max0, base::next_index1, Max1, base::next_index2, Max2, Assign>  next_t;

    template <typename MatrixA, typename MatrixB, typename MatrixC>
    static inline void apply(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	Assign::update(C(base::index1, base::index2), A(base::index1, base::index0) * B(base::index0, base::index2));
	next_t::apply(A, B, C);
    }   
};

template <std::size_t Max0, std::size_t Max1, std::size_t Max2, typename Assign>
struct fully_unroll_dmat_dmat_mult_block<Max0, Max0, Max1, Max1, Max2, Max2, Assign>
  : public meta_math::loop3<Max0, Max0, Max1, Max1, Max2, Max2>
{
    typedef meta_math::loop3<Max0, Max0, Max1, Max1, Max2, Max2>                              base;

    template <typename MatrixA, typename MatrixB, typename MatrixC>
    static inline void apply(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {
	Assign::update(C(base::index1, base::index2), A(base::index1, base::index0) * B(base::index0, base::index2));
    }   
};

template <typename Assign= assign::assign_sum, 
	  typename Backup= gen_dmat_dmat_mult_t<Assign> >
struct fully_unroll_fixed_size_dmat_dmat_mult_t
{
    // if C is empty just do nothing
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    typename boost::enable_if_c<static_size<MatrixC>::value == 0>::type
    operator()(MatrixA const&, MatrixB const&, MatrixC&) {}

    // just initialize
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    typename boost::enable_if_c<static_num_cols<MatrixA>::value == 0 && static_size<MatrixC>::value != 0>::type
    operator()(MatrixA const&, MatrixB const&, MatrixC& C) { if (Assign::init_to_zero) set_to_zero(C); }

    struct noop
    {
	template <typename MatrixA, typename MatrixB, typename MatrixC>
	static inline void apply(MatrixA const&, MatrixB const&, MatrixC&) {}
    };

    // really compute
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    typename boost::enable_if_c<static_size<MatrixA>::value != 0 && static_size<MatrixC>::value != 0>::type
    operator()(MatrixA const& A, MatrixB const& B, MatrixC& C)
    {	
	vampir_trace<1002> tracer;
	typedef typename static_num_rows<MatrixC>::type size_type;
	static const size_type rows_c= static_num_rows<MatrixC>::value, cols_c= static_num_cols<MatrixC>::value, 
	  	               cols_a= static_num_cols<MatrixA>::value;
	// corresponds to C= A[all][0] * B[0][all];
	fully_unroll_dmat_dmat_mult_init_block<1, rows_c, 1, cols_c, Assign>::apply(A, B, C); 

	// corresponds to C+= A[all][1:] * B[1:][all]; if necessary
	typedef fully_unroll_dmat_dmat_mult_block<2, cols_a, 1, rows_c, 1, cols_c, Assign>  f2;
	typedef typename boost::mpl::if_c<(cols_a > 1), f2, noop>::type                     f3;
	f3::apply(A, B, C);
    }
};    


} // namespace mtl

#endif // MTL_DMAT_DMAT_MULT_INCLUDE

// Include plattform specific implementations
#include <boost/numeric/mtl/operation/opteron/matrix_mult.hpp>

