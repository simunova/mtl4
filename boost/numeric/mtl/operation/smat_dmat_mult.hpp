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

#ifndef MTL_SMAT_DMAT_MULT_INCLUDE
#define MTL_SMAT_DMAT_MULT_INCLUDE

#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/flatcat.hpp>
#include <boost/numeric/meta_math/loop1.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace functor {

template <typename Assign= assign::assign_sum,
	  typename Backup= no_op>     // To allow 2nd parameter, is ignored
struct gen_smat_dmat_mult
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& a, MatrixB const& b, MatrixC& c)
    {
	vampir_trace<4018> tracer;
	apply(a, b, c, typename OrientedCollection<MatrixA>::orientation());
    }

private:
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& a, MatrixB const& b, MatrixC& c, tag::row_major)
    {
	using namespace tag;
	using traits::range_generator;  
        typedef typename range_generator<row, MatrixA>::type       a_cur_type;             
        typedef typename range_generator<row, MatrixC>::type       c_cur_type;             
	typedef typename range_generator<col, MatrixB>::type       b_cur_type;    
         
        typedef typename range_generator<nz, a_cur_type>::type     a_icur_type;            
        typedef typename range_generator<all, b_cur_type>::type    b_icur_type;          
        typedef typename range_generator<iter::all, c_cur_type>::type    c_icur_type;            

	typename traits::col<MatrixA>::type             col_a(a); 
	typename traits::const_value<MatrixA>::type     value_a(a); 
	typename traits::const_value<MatrixB>::type     value_b(b); 

	if (Assign::init_to_zero) set_to_zero(c);

	a_cur_type ac= begin<row>(a), aend= end<row>(a);
	for (c_cur_type cc= begin<row>(c); ac != aend; ++ac, ++cc) {

	    b_cur_type bc= begin<col>(b), bend= end<col>(b);
	    for (c_icur_type cic= begin<iter::all>(cc); bc != bend; ++bc, ++cic) { 
		    
		typename MatrixC::value_type c_tmp(*cic);
		for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) {

		    typename Collection<MatrixA>::size_type     ca= col_a(*aic);   // column of non-zero

		    b_icur_type bic= begin<all>(bc);
		    bic+= ca;
		    Assign::update(c_tmp, value_a(*aic) * value_b(*bic));
		}
		*cic= c_tmp;
	    }
	}
    }


    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& a, MatrixB const& b, MatrixC& c, tag::col_major)
    {
	using namespace tag;
	using traits::range_generator;  
        typedef typename range_generator<col, MatrixA>::type       a_cur_type;             
        typedef typename range_generator<nz, a_cur_type>::type     a_icur_type;            

	typename traits::row<MatrixA>::type             row_a(a); 
	typename traits::const_value<MatrixA>::type     value_a(a); 

	if (Assign::init_to_zero) set_to_zero(c);

	unsigned rb= 0; // traverse all rows of b
	for (a_cur_type ac= begin<col>(a), aend= end<col>(a); ac != aend; ++ac, ++rb)
	    for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) {
		typename Collection<MatrixA>::size_type     ra= row_a(*aic);   // row in A and C
		typename Collection<MatrixA>::value_type    va= value_a(*aic); // value of non-zero

		for (unsigned cb= 0; cb < num_cols(b); ++cb) // column in B and C
		    Assign::update(c(ra, cb), va * b(rb, cb));
	    }
    }
};


// =======================
// Unrolled 
// required has_2D_layout
// =======================

// Define defaults if not yet given as Compiler flag
#ifndef MTL_SMAT_DMAT_MULT_TILING1
#  define MTL_SMAT_DMAT_MULT_TILING1 8
#endif

template <unsigned long Index0, unsigned long Max0, typename Assign>
struct gen_tiling_smat_dmat_mult_block
    : public meta_math::loop1<Index0, Max0>
{
    typedef meta_math::loop1<Index0, Max0>                                    base;
    typedef gen_tiling_smat_dmat_mult_block<base::next_index0, Max0, Assign>  next_t;

    template <typename Value, typename ValueA, typename ValueB, typename Size>
    static inline void apply(Value& tmp00, Value& tmp01, Value& tmp02, Value& tmp03, Value& tmp04, 
			     Value& tmp05, Value& tmp06, Value& tmp07, Value& tmp08, Value& tmp09, 
			     Value& tmp10, Value& tmp11, Value& tmp12, Value& tmp13, Value& tmp14, Value& tmp15, 
			     const ValueA& va, ValueB *begin_b, const Size& bci)
    {
	tmp00+= va * *(begin_b + base::index0 * bci);
	next_t::apply(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		      tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
		      va, begin_b, bci); 
    }

    template <typename Value, typename MatrixC, typename SizeC>
    static inline void update(Value& tmp00, Value& tmp01, Value& tmp02, Value& tmp03, Value& tmp04, 
			      Value& tmp05, Value& tmp06, Value& tmp07, Value& tmp08, Value& tmp09, 
			      Value& tmp10, Value& tmp11, Value& tmp12, Value& tmp13, Value& tmp14, Value& tmp15,
			      MatrixC& c, SizeC i, SizeC k)
    {
	Assign::update(c(i, k + base::index0), tmp00);
	next_t::update(tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
		       tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp00, 
		       c, i, k);
    }
};
	

template <unsigned long Max0, typename Assign>
struct gen_tiling_smat_dmat_mult_block<Max0, Max0, Assign>
    : public meta_math::loop1<Max0, Max0>
{
    typedef meta_math::loop1<Max0, Max0>                                    base;

    template <typename Value, typename ValueA, typename ValueB, typename Size>
    static inline void apply(Value& tmp00, Value&, Value&, Value&, Value&, 
			     Value&, Value&, Value&, Value&, Value&, 
			     Value&, Value&, Value&, Value&, Value&, Value&, 
			     const ValueA& va, ValueB *begin_b, const Size& bci)
    {
	tmp00+= va * *(begin_b + base::index0 * bci);
    }

    template <typename Value, typename MatrixC, typename SizeC>
    static inline void update(Value& tmp00, Value&, Value&, Value&, Value&, 
			      Value&, Value&, Value&, Value&, Value&, 
			      Value&, Value&, Value&, Value&, Value&, Value&,
			      MatrixC& c, SizeC i, SizeC k)
    {
	Assign::update(c(i, k + base::index0), tmp00);
    }
};


template <unsigned long Tiling1= MTL_SMAT_DMAT_MULT_TILING1,
	  typename Assign= assign::assign_sum,
	  typename Backup= gen_smat_dmat_mult<Assign> >
struct gen_tiling_smat_dmat_mult
{
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void operator()(MatrixA const& a, MatrixB const& b, MatrixC& c)
    {
    vampir_trace<4019> tracer;
	apply(a, b, c, traits::layout_flatcat<MatrixC>());
    }

private:
    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& a, MatrixB const& b, MatrixC& c, tag::universe)
    {
	Backup()(a, b, c);
    }

    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply(MatrixA const& a, MatrixB const& b, MatrixC& c, tag::flat<tag::has_2D_layout>)
    {
	apply2(a, b, c, typename OrientedCollection<MatrixA>::orientation());
    }

    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply2(MatrixA const& a, MatrixB const& b, MatrixC& c, tag::col_major)
    {
	// may be I'll write an optimized version later
	Backup()(a, b, c);
    }


    template <typename MatrixA, typename MatrixB, typename MatrixC>
    void apply2(MatrixA const& a, MatrixB const& b, MatrixC& c, tag::row_major)
    {
	using namespace tag;
	using traits::range_generator;  

	typedef gen_tiling_smat_dmat_mult_block<1, Tiling1, Assign>  block;
	typedef typename Collection<MatrixA>::size_type              size_type;
	typedef typename Collection<MatrixC>::value_type             value_type;
	const value_type z= math::zero(c[0][0]);    // if this are matrices we need their size

        typedef typename range_generator<row, MatrixA>::type         a_cur_type;             
        typedef typename range_generator<nz, a_cur_type>::type       a_icur_type;            

	typename traits::col<MatrixA>::type             col_a(a); 
	typename traits::const_value<MatrixA>::type     value_a(a); 

	if (Assign::init_to_zero) set_to_zero(c);

	size_type i_max= num_cols(b), i_block= Tiling1 * (i_max / Tiling1);
	size_t bci= i_max > 1 ? &b(0, 1) - &b(0, 0) : 1; // offset of incrementing B's column if more than 1 column

	size_type rc= 0; // start in row 0
	for (a_cur_type ac= begin<row>(a), aend= end<row>(a); ac != aend; ++ac, ++rc) {

	    for (size_type i= 0; i < i_block; i+= Tiling1) {
	    
		value_type tmp00= z, tmp01= z, tmp02= z, tmp03= z, tmp04= z,
                           tmp05= z, tmp06= z, tmp07= z, tmp08= z, tmp09= z,
 		           tmp10= z, tmp11= z, tmp12= z, tmp13= z, tmp14= z, tmp15= z;

		for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) {
		    typename Collection<MatrixA>::size_type     ca= col_a(*aic);   // column of non-zero
		    typename Collection<MatrixA>::value_type    va= value_a(*aic); // value of non-zero

		    // Element in first vector in block to be multiplied with va; rb==ca
		    const typename MatrixB::value_type *begin_b= &b(ca, i); 
		    block::apply(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
				 tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
				 va, begin_b, bci); 
		}
		block::update(tmp00, tmp01, tmp02, tmp03, tmp04, tmp05, tmp06, tmp07, tmp08, tmp09, 
			      tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, 
			      c, rc, i);
	    }	

	    for (size_type i= i_block; i < i_max; i++) {
		value_type tmp00= z;
		for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) {
		    typename Collection<MatrixA>::size_type     ca= col_a(aic);   // column of non-zero
		    tmp00+= value_a(*aic) * b(ca, i);
		}
		Assign::update(c(rc, i), tmp00);
	    }
	}
    }
};




}} // namespace mtl::functor

#endif // MTL_SMAT_DMAT_MULT_INCLUDE
