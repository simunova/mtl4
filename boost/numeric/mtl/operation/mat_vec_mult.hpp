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

#ifndef MTL_MAT_VEC_MULT_INCLUDE
#define MTL_MAT_VEC_MULT_INCLUDE

#include <cassert>
// #include <iostream>
#include <boost/mpl/bool.hpp>

#include <boost/numeric/mtl/config.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/multi_tmp.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/utility/omp_size_type.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/meta_math/loop.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace mat {

namespace impl {

    template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1, typename Assign>
    struct fully_unroll_mat_cvec_mult
      : public meta_math::loop2<Index0, Max0, Index1, Max1>
    {
	typedef meta_math::loop2<Index0, Max0, Index1, Max1>                              base;
	typedef fully_unroll_mat_cvec_mult<base::next_index0, Max0, base::next_index1, Max1, Assign>  next_t;

	template <typename Matrix, typename VectorIn, typename VectorOut>
	static inline void apply(const Matrix& A, const VectorIn& v, VectorOut& w)
	{
	    Assign::update(w[base::index0], A[base::index0][base::index1] * v[base::index1]);
	    next_t::apply(A, v, w);
	}   
    };

    // need specialization here for not going back to column 0 but column 1
    template <std::size_t Index0, std::size_t Max0, std::size_t Max1, typename Assign>
    struct fully_unroll_mat_cvec_mult<Index0, Max0, Max1, Max1, Assign>
      : public meta_math::loop2<Index0, Max0, Max1, Max1>
    {
	typedef meta_math::loop2<Index0, Max0, Max1, Max1>                              base;
	typedef fully_unroll_mat_cvec_mult<base::next_index0, Max0, 2, Max1, Assign>  next_t;

	template <typename Matrix, typename VectorIn, typename VectorOut>
	static inline void apply(const Matrix& A, const VectorIn& v, VectorOut& w)
	{
	    Assign::update(w[base::index0], A[base::index0][base::index1] * v[base::index1]);
	    next_t::apply(A, v, w);
	}   
    };

    template <std::size_t Max0, std::size_t Max1, typename Assign>
    struct fully_unroll_mat_cvec_mult<Max0, Max0, Max1, Max1, Assign>
      : public meta_math::loop2<Max0, Max0, Max1, Max1>
    {
	typedef meta_math::loop2<Max0, Max0, Max1, Max1>                              base;

	template <typename Matrix, typename VectorIn, typename VectorOut>
	static inline void apply(const Matrix& A, const VectorIn& v, VectorOut& w)
	{
	    Assign::update(w[base::index0], A[base::index0][base::index1] * v[base::index1]);
	}   
    };

    struct noop
    {
	template <typename Matrix, typename VectorIn, typename VectorOut>
	static inline void apply(const Matrix&, const VectorIn&, VectorOut&) {}
    };
} // impl

// Dense matrix vector multiplication with run-time matrix size
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void dense_mat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, boost::mpl::true_)
{
    vampir_trace<3017> tracer;
    typedef typename static_num_rows<Matrix>::type size_type;
    static const size_type rows_a= static_num_rows<Matrix>::value, cols_a= static_num_cols<Matrix>::value;

    assert(rows_a > 0 && cols_a > 0);
    // w= A[all][0] * v[0];  N.B.: 1D is unrolled by the compiler faster (at least on gcc)
    for (size_type i= 0; i < rows_a; i++) 
	Assign::first_update(w[i], A[i][0] * v[0]);
	
    // corresponds to w+= A[all][1:] * v[1:]; if necessary
    typedef impl::fully_unroll_mat_cvec_mult<1, rows_a, 2, cols_a, Assign>  f2;
    typedef typename boost::mpl::if_c<(cols_a > 1), f2, impl::noop>::type   f3;
    f3::apply(A, v, w);
}
	
template <unsigned Index, unsigned Size>
struct init_ptrs
{
    template <typename Matrix, typename Ptrs>
    inline static void apply(const Matrix& A, Ptrs& ptrs)
    {
	ptrs.value= &A[Index][0];
	init_ptrs<Index+1, Size>::apply(A, ptrs.sub);
    }
};

template <unsigned Size>
struct init_ptrs<Size, Size>
{
    template <typename Matrix, typename Ptrs>
    inline static void apply(const Matrix&, Ptrs&) {}
};

template <unsigned Index, unsigned Size>
struct square_cvec_mult_rows
{
    template <typename Tmps, typename Ptrs, typename ValueIn>
    inline static void compute(Tmps& tmps, Ptrs& ptrs, ValueIn vi)
    {
	tmps.value+= *ptrs.value++ * vi;
	square_cvec_mult_rows<Index+1, Size>::compute(tmps.sub, ptrs.sub, vi);
    }

    template <typename As, typename VectorOut, typename Tmps>
    inline static void update(As, VectorOut& w, const Tmps& tmps)
    {
	As::first_update(w[Index], tmps.value);
	square_cvec_mult_rows<Index+1, Size>::update(As(), w, tmps.sub);
    }
};

template <unsigned Size>
struct square_cvec_mult_rows<Size, Size>
{
    template <typename Tmps, typename Ptrs, typename ValueIn>
    inline static void compute(Tmps&, Ptrs&, ValueIn) {}

    template <typename As, typename VectorOut, typename Tmps>
    inline static void update(As, VectorOut&, const Tmps&) {}
};

template <unsigned Index, unsigned Size>
struct square_cvec_mult_cols
{    
    template <typename Tmps, typename Ptrs, typename VPtr>
    inline static void compute(Tmps& tmps, Ptrs& ptrs, VPtr vp)
    {
	square_cvec_mult_rows<0, Size>::compute(tmps, ptrs, *vp);
	square_cvec_mult_cols<Index+1, Size>::compute(tmps, ptrs, ++vp);
    }
};

template <unsigned Size>
struct square_cvec_mult_cols<Size, Size>
{    
    template <typename Tmps, typename Ptrs, typename VPtr>
    inline static void compute(Tmps&, Ptrs&, VPtr) {}
};


// Dense matrix vector multiplication with run-time matrix size
template <unsigned Size, typename MValue, typename MPara, typename ValueIn, typename ParaIn, 
	  typename VectorOut, typename Assign>
inline void square_cvec_mult(const dense2D<MValue, MPara>& A, const mtl::vec::dense_vector<ValueIn, ParaIn>& v, VectorOut& w, Assign)
{
    // vampir_trace<3067> tracer;
    MTL_STATIC_ASSERT((mtl::traits::is_row_major<MPara>::value), "Only row-major matrices supported in this function.");

    typedef typename Collection<VectorOut>::value_type value_type;    
    multi_tmp<Size, value_type> tmps(math::zero(w[0]));

    multi_tmp<Size, const MValue*>    ptrs;
    init_ptrs<0, Size>::apply(A, ptrs);

    const ValueIn* vp= &v[0];

    square_cvec_mult_cols<0, Size>::compute(tmps, ptrs, vp); // outer loop over columns
    square_cvec_mult_rows<0, Size>::update(Assign(), w, tmps);         // update rows
}


// Dense matrix vector multiplication with run-time matrix size
template <typename MValue, typename MPara, typename ValueIn, typename ParaIn, typename VectorOut, typename Assign>
typename boost::enable_if<mtl::traits::is_row_major<MPara> >::type
inline dense_mat_cvec_mult(const dense2D<MValue, MPara>& A, const mtl::vec::dense_vector<ValueIn, ParaIn>& v, VectorOut& w, Assign, boost::mpl::false_)
{
    // vampir_trace<3066> tracer;

    using math::zero; 
    if (mtl::size(w) == 0) return;

    typedef typename Collection<VectorOut>::value_type value_type;
    // typedef ValueIn                                    value_in_type;
    typedef typename MPara::size_type                  size_type;

    const size_type  nr= num_rows(A), nc= num_cols(A);
    if (nr == nc && nr <= 8) 
	switch (nr) {
	  case 1: Assign::first_update(w[0], A[0][0] * v[0]); return;
	  case 2: square_cvec_mult<2>(A, v, w, Assign()); return;
	  case 3: square_cvec_mult<3>(A, v, w, Assign()); return;
	  case 4: square_cvec_mult<4>(A, v, w, Assign()); return;
	  case 5: square_cvec_mult<5>(A, v, w, Assign()); return;
	  case 6: square_cvec_mult<6>(A, v, w, Assign()); return;
	  case 7: square_cvec_mult<7>(A, v, w, Assign()); return;
	  case 8: square_cvec_mult<8>(A, v, w, Assign()); return;
	}


    const size_type  nrb= nr / 4 * 4;
    const value_type z(math::zero(w[0]));

    for (size_type i= 0; i < nrb; i+= 4) {
	value_type      tmp0(z), tmp1(z), tmp2(z), tmp3(z);
	const MValue *p0= &A[i][0], *pe= p0 + nc, *p1= &A[i+1][0], *p2= &A[i+2][0], *p3= &A[i+3][0];
	const ValueIn* vp= &v[0];
	for (; p0 != pe; ) {
	    const ValueIn vj= *vp++;
	    tmp0+= *p0++ * vj;
	    tmp1+= *p1++ * vj;
	    tmp2+= *p2++ * vj;
	    tmp3+= *p3++ * vj;
	}
	Assign::first_update(w[i], tmp0);
	Assign::first_update(w[i+1], tmp1);
	Assign::first_update(w[i+2], tmp2);
	Assign::first_update(w[i+3], tmp3);
    }

    for (size_type i= nrb; i < nr; i++) {
	value_type tmp= z;
	const ValueIn* vp= &v[0];
	for (const MValue *p0= &A[i][0], *pe= p0 + nc; p0 != pe; ) 
	    tmp+= *p0++ * *vp++;
	Assign::first_update(w[i], tmp);
    }
}


// Dense matrix vector multiplication with run-time matrix size
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void dense_mat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, boost::mpl::false_)
{
    vampir_trace<3018> tracer;
    // Naive implementation, will be moved to a functor and complemented with more efficient ones

    using math::zero; 
    if (mtl::size(w) == 0) return;
    // std::cout << "Bin in richtiger Funktion\n";

    // if (Assign::init_to_zero) set_to_zero(w); // replace update with first_update insteda

    typedef typename Collection<VectorOut>::value_type value_type;
    typedef typename Collection<VectorIn>::value_type  value_in_type;
    typedef typename Collection<Matrix>::size_type     size_type;

    const value_type z(math::zero(w[0]));
    const size_type nr= num_rows(A), nrb= nr / 4 * 4, nc= num_cols(A);

    for (size_type i= 0; i < nrb; i+= 4) {
	value_type      tmp0(z), tmp1(z), tmp2(z), tmp3(z);
	for (size_type j= 0; j < nc; j++) {
	    const value_in_type vj= v[j];
	    tmp0+= A[i][j] * vj;
	    tmp1+= A[i+1][j] * vj;
	    tmp2+= A[i+2][j] * vj;
	    tmp3+= A[i+3][j] * vj;
	}
	Assign::first_update(w[i], tmp0);
	Assign::first_update(w[i+1], tmp1);
	Assign::first_update(w[i+2], tmp2);
	Assign::first_update(w[i+3], tmp3);
    }

    for (size_type i= nrb; i < nr; i++) {
	value_type tmp= z;
	for (size_type j= 0; j < nc; j++) 
	    tmp+= A[i][j] * v[j];
	Assign::first_update(w[i], tmp);
    }
}

// Dense matrix vector multiplication
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void mat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::dense>)
{
# ifdef MTL_NOT_UNROLL_FSIZE_MAT_VEC_MULT
    boost::mpl::false_        selector;
# else
	mtl::traits::is_static<Matrix> selector;
# endif
    dense_mat_cvec_mult(A, v, w, Assign(), selector);
}

// Element structure vector multiplication
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void mat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::element_structure>)
{
    vampir_trace<3048> tracer;
    if (mtl::size(w) == 0) return;

    typedef typename Collection<VectorOut>::value_type value_type;
    typedef typename Collection<VectorIn>::value_type  value_in_type;

    value_in_type varray[1024];
    value_type    warray[1024];

    if (Assign::init_to_zero) set_to_zero(w);
    for(int elmi= 0; elmi < A.m_total_elements; elmi++){
	const typename Matrix::element_type& elementi= A.m_elements[elmi];
	const typename Matrix::element_type::index_type& indices= elementi.get_indices();
	std::size_t n= size(indices);

	if (n <= 1024) {
	    VectorIn vtmp(n, varray);
	    for (unsigned int i= 0; i < n; i++)
		vtmp[i]= v[indices[i]];
	    VectorOut wtmp(n, warray);
	    wtmp= elementi.get_values() * vtmp;
	    for (unsigned int i= 0; i < n; i++)
		Assign::update(w[indices[i]], wtmp[i]);
	} else {
	    VectorIn vtmp(n);
	    for (unsigned int i= 0; i < n; i++)
		vtmp[i]= v[indices[i]];
	    VectorOut wtmp(elementi.get_values() * vtmp);
	    for (unsigned int i= 0; i < n; i++)
		Assign::update(w[indices[i]], wtmp[i]);
	}
    }
}

// Multi-vector vector multiplication (tag::multi_vector is derived from dense)
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void mat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::multi_vector>)
{
    vampir_trace<3019> tracer;
    if (Assign::init_to_zero) set_to_zero(w);
    for (unsigned i= 0; i < num_cols(A); i++)
	Assign::update(w, A.vector(i) * v[i]);
}

// Transposed multi-vector vector multiplication (tag::transposed_multi_vector is derived from dense)
template <typename TransposedMatrix, typename VectorIn, typename VectorOut, typename Assign>
inline void mat_cvec_mult(const TransposedMatrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::transposed_multi_vector>)
{
    vampir_trace<3020> tracer;
    typename TransposedMatrix::const_ref_type B= A.ref; // Referred matrix

    if (Assign::init_to_zero) set_to_zero(w);
    for (unsigned i= 0; i < num_cols(B); i++)
	Assign::update(w[i], dot_real(B.vector(i), v));
}

// Hermitian multi-vector vector multiplication (tag::hermitian_multi_vector is derived from dense)
template <typename HermitianMatrix, typename VectorIn, typename VectorOut, typename Assign>
inline void mat_cvec_mult(const HermitianMatrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::hermitian_multi_vector>)
{
    vampir_trace<3021> tracer;
    typename HermitianMatrix::const_ref_type B= A.const_ref(); // Referred matrix

    if (Assign::init_to_zero) set_to_zero(w);
    for (unsigned i= 0; i < num_cols(B); i++)
	Assign::update(w[i], dot(B.vector(i), v));
}



// Sparse row-major matrix vector multiplication
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void smat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::row_major)
{
    vampir_trace<3022> tracer;
    using namespace tag; 
    using mtl::traits::range_generator;  
    using math::zero;
    using mtl::vec::set_to_zero;

    typedef typename range_generator<row, Matrix>::type       a_cur_type;    
    typedef typename range_generator<nz, a_cur_type>::type    a_icur_type;            
    typename mtl::traits::col<Matrix>::type                   col_a(A); 
    typename mtl::traits::const_value<Matrix>::type           value_a(A); 

    if (Assign::init_to_zero) set_to_zero(w);

    typedef typename Collection<VectorOut>::value_type        value_type;
    a_cur_type ac= begin<row>(A), aend= end<row>(A);
    for (unsigned i= 0; ac != aend; ++ac, ++i) {
	value_type tmp= zero(w[i]);
	for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) 
	    tmp+= value_a(*aic) * v[col_a(*aic)];	
	Assign::update(w[i], tmp);
    }
}

// Row-major compressed2D with very few entries (i.e. Very Sparse MATrix) times vector
template <typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
inline void vsmat_cvec_mult(const compressed2D<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign, tag::row_major)
{
    vampir_trace<3064> tracer;
    using math::zero;

    typedef compressed2D<MValue, MPara>                       Matrix;
    typedef typename Collection<Matrix>::size_type            size_type; 
    typedef typename Collection<VectorOut>::value_type        value_type;

    if (mtl::size(w) == 0) return;
    const value_type z(math::zero(w[0]));

    // std::cout << "very sparse: nnz = " << A.nnz() << ", num_rows = " << num_rows(A) << '\n';

    size_type nr= num_rows(A);
    for (size_type i1= 0, i2= std::min<size_type>(1024, nr); i1 < i2; i1= i2, i2= std::min<size_type>(i2 + 1024, nr)) {
	// std::cout << "range = " << i1 << " .. " << i2 << "\n";
	if (A.ref_major()[i1] < A.ref_major()[i2])
	    for (size_type i= i1; i < i2; ++i) {
		const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1];
		value_type      tmp0(z);
		for (size_type j0= cj0; j0 != cj1; ++j0)
		    tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
		Assign::first_update(w[i], tmp0);
	    }
    }
}

#ifdef MTL_CRS_CVEC_MULT_TUNING
template <unsigned Index, unsigned BSize, typename SizeType>
struct crs_cvec_mult_block
{
    template <typename Matrix, typename VectorIn, typename CBlock, typename TBlock>
    void operator()(const Matrix& A, const VectorIn& v, const CBlock& cj, TBlock& tmp) const
    {
	for (SizeType j= cj.value; j != cj.sub.value; ++j) // cj is one index larger
	    tmp.value+= A.data[j] * v[A.ref_minor()[j]];
	sub(A, v, cj.sub, tmp.sub);
    }

    template <typename VectorOut, typename TBlock, typename Assign>
    void first_update(VectorOut& w, SizeType i, const TBlock& tmp, Assign as) const
    { 
	Assign::first_update(w[i + Index], tmp.value);
	sub.first_update(w, i, tmp.sub, as);
    }
    
    crs_cvec_mult_block<Index+1, BSize, SizeType> sub;
};


template <unsigned BSize, typename SizeType>
struct crs_cvec_mult_block<BSize, BSize, SizeType>
{
    template <typename Matrix, typename VectorIn, typename CBlock, typename TBlock>
    void operator()(const Matrix& A, const VectorIn& v, const CBlock& cj, TBlock& tmp) const
    {
	for (SizeType j= cj.value; j != cj.sub.value; ++j)// cj is one index larger
	    tmp.value+= A.data[j] * v[A.ref_minor()[j]];
    }

    template <typename VectorOut, typename TBlock, typename Assign>
    void first_update(VectorOut& w, SizeType i, const TBlock& tmp, Assign) const
    { 
	Assign::first_update(w[i + BSize], tmp.value);
    }
};


// Row-major compressed2D vector multiplication
template <unsigned BSize, typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
inline void smat_cvec_mult(const compressed2D<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign as, tag::row_major)
{
    vampir_trace<3049> tracer;
    using math::zero;

    if (A.nnz() < num_rows(A)) {
	vsmat_cvec_mult(A, v, w, as, tag::row_major());
	return;
    }

    typedef compressed2D<MValue, MPara>                       Matrix;
    typedef typename Collection<VectorOut>::value_type        value_type;
    typedef typename mtl::traits::omp_size_type<typename Collection<Matrix>::size_type>::type size_type;

    if (size(w) == 0) return;
    const value_type z(math::zero(w[0]));

    size_type nr= num_rows(A), nrb= nr / BSize * BSize;

    #ifdef MTL_WITH_OPENMP
    #   pragma omp parallel
    #endif
    {
    	#ifdef MTL_WITH_OPENMP
	    vampir_trace<8004> tracer;
    	#   pragma omp for
    	#endif
	for (size_type i= 0; i < nrb; i+= BSize) {
	    multi_constant_from_array<0, BSize+1, size_type> cj(A.ref_major(), i);
	    multi_tmp<BSize, value_type>                     tmp(z);
	    crs_cvec_mult_block<0, BSize-1, size_type>       block;

	    block(A, v, cj, tmp);
	    block.first_update(w, i, tmp, as);
	}
    }

    for (size_type i= nrb; i < nr; ++i) {
	const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1];
	value_type      tmp0(z);
	for (size_type j0= cj0; j0 != cj1; ++j0)
	    tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
	Assign::first_update(w[i], tmp0);
    }
}

template <typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
typename mtl::traits::enable_if_scalar<typename Collection<VectorOut>::value_type>::type
inline smat_cvec_mult(const compressed2D<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign, tag::row_major)
{
    smat_cvec_mult<crs_cvec_mult_block_size>(A, v, w, Assign(), tag::row_major());
}
#endif


#if !defined(MTL_CRS_CVEC_MULT_NO_ACCEL) && !defined(MTL_CRS_CVEC_MULT_TUNING)

template <typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
typename mtl::traits::enable_if_scalar<typename Collection<VectorOut>::value_type>::type
inline adapt_crs_cvec_mult(const compressed2D<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign)
{
    vampir_trace<3065> tracer;
    using math::zero;
    assert(!Assign::init_to_zero);

    typedef compressed2D<MValue, MPara>                       Matrix;
    typedef typename Collection<Matrix>::size_type            size_type; 
    typedef typename Collection<VectorOut>::value_type        value_type;

    const value_type z(math::zero(w[0]));
    size_type nr= num_rows(A), nrb= nr / 4 * 4, nrb2= nr / 64 * 64;

    for (size_type i1= 0; i1 < nrb2; i1+= 64) 
	if (A.ref_major()[i1] != A.ref_major()[i1 + 64])
	    for (size_type i2= i1, i2e= i1+64; i2 < i2e; i2+= 16)
		if (A.ref_major()[i2] != A.ref_major()[i2 + 16])
		    for (size_type i= i2, i3e= i2+16; i < i3e; i+= 4) 
			if (A.ref_major()[i] != A.ref_major()[i + 4]) {
			    const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1], cj2= A.ref_major()[i+2], 
				cj3= A.ref_major()[i+3], cj4= A.ref_major()[i+4];
			    value_type      tmp0(z), tmp1(z), tmp2(z), tmp3(z);
			    for (size_type j0= cj0; j0 != cj1; ++j0)
				tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
			    for (size_type j1= cj1; j1 != cj2; ++j1)
				tmp1+= A.data[j1] * v[A.ref_minor()[j1]];
			    for (size_type j2= cj2; j2 != cj3; ++j2)
				tmp2+= A.data[j2] * v[A.ref_minor()[j2]];
			    for (size_type j3= cj3; j3 != cj4; ++j3)
				tmp3+= A.data[j3] * v[A.ref_minor()[j3]];

			    Assign::first_update(w[i], tmp0);
			    Assign::first_update(w[i+1], tmp1);
			    Assign::first_update(w[i+2], tmp2);
			    Assign::first_update(w[i+3], tmp3);
			}

    for (size_type i= nrb2; i < nrb; i+= 4) 
	if (A.ref_major()[i] != A.ref_major()[i + 4])  {
	    const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1], cj2= A.ref_major()[i+2], 
		cj3= A.ref_major()[i+3], cj4= A.ref_major()[i+4];
	    value_type      tmp0(z), tmp1(z), tmp2(z), tmp3(z);
	    for (size_type j0= cj0; j0 != cj1; ++j0)
		tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
	    for (size_type j1= cj1; j1 != cj2; ++j1)
		tmp1+= A.data[j1] * v[A.ref_minor()[j1]];
	    for (size_type j2= cj2; j2 != cj3; ++j2)
		tmp2+= A.data[j2] * v[A.ref_minor()[j2]];
	    for (size_type j3= cj3; j3 != cj4; ++j3)
		tmp3+= A.data[j3] * v[A.ref_minor()[j3]];

	    Assign::first_update(w[i], tmp0);
	    Assign::first_update(w[i+1], tmp1);
	    Assign::first_update(w[i+2], tmp2);
	    Assign::first_update(w[i+3], tmp3);
	}

    for (size_type i= nrb; i < nr; ++i) {
	const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1];
	value_type      tmp0(z);
	for (size_type j0= cj0; j0 != cj1; ++j0)
	    tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
	Assign::first_update(w[i], tmp0);
    }
}

// Row-major compressed2D vector multiplication
template <typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
typename mtl::traits::enable_if_scalar<typename Collection<VectorOut>::value_type>::type
inline smat_cvec_mult(const compressed2D<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign as, tag::row_major)
{
    vampir_trace<3049> tracer;
    // vampir_trace<5056> tttracer;
    using math::zero;

    if (A.nnz() < num_rows(A) && !as.init_to_zero) {
	vsmat_cvec_mult(A, v, w, as, tag::row_major());
	return;
    }

    typedef compressed2D<MValue, MPara>                       Matrix;
    typedef typename Collection<VectorOut>::value_type        value_type;
    typedef typename mtl::traits::omp_size_type<typename Collection<Matrix>::size_type>::type size_type;

    if (mtl::size(w) == 0) return;
    const value_type z(math::zero(w[0]));

    size_type nr= num_rows(A), nrb= nr / 4 * 4;
    if (nr > 10) {
	size_type nh= nr / 2, nq= nr / 4, nt= nr - nq;
	if (!as.init_to_zero &&
	    (A.ref_major()[1] == A.ref_major()[0] 
	     || A.ref_major()[nq] == A.ref_major()[nq+1]
	     || A.ref_major()[nh] == A.ref_major()[nh+1]
	     || A.ref_major()[nt] == A.ref_major()[nt+1]
	     || A.ref_major()[nr-1] == A.ref_major()[nr])) {
	    adapt_crs_cvec_mult(A, v, w, as);
	    return;
	}
    }

    #ifdef MTL_WITH_OPENMP
    #   pragma omp parallel
    #endif
    {
    	#ifdef MTL_WITH_OPENMP
	    vampir_trace<8004> tracer;
    	#   pragma omp for
    	#endif
	    for (size_type i= 0; i < nrb; i+= 4) {
		const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1], cj2= A.ref_major()[i+2], 
		    cj3= A.ref_major()[i+3], cj4= A.ref_major()[i+4];
		value_type      tmp0(z), tmp1(z), tmp2(z), tmp3(z);
		for (size_type j0= cj0; j0 != cj1; ++j0)
		    tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
		for (size_type j1= cj1; j1 != cj2; ++j1)
		    tmp1+= A.data[j1] * v[A.ref_minor()[j1]];
		for (size_type j2= cj2; j2 != cj3; ++j2)
		    tmp2+= A.data[j2] * v[A.ref_minor()[j2]];
		for (size_type j3= cj3; j3 != cj4; ++j3)
		    tmp3+= A.data[j3] * v[A.ref_minor()[j3]];

		Assign::first_update(w[i], tmp0);
		Assign::first_update(w[i+1], tmp1);
		Assign::first_update(w[i+2], tmp2);
		Assign::first_update(w[i+3], tmp3);
	    }
    }

    for (size_type i= nrb; i < nr; ++i) {
	const size_type cj0= A.ref_major()[i], cj1= A.ref_major()[i+1];
	value_type      tmp0(z);
	for (size_type j0= cj0; j0 != cj1; ++j0)
	    tmp0+= A.data[j0] * v[A.ref_minor()[j0]];
	Assign::first_update(w[i], tmp0);
    }
}
#endif

// Row-major ell_matrix vector multiplication
template <typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
typename mtl::traits::enable_if_scalar<typename Collection<VectorOut>::value_type>::type
inline smat_cvec_mult(const ell_matrix<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign, tag::row_major)
{
    typedef typename MPara::size_type size_type;

    const size_type stride= A.stride(), slots= A.slots();
    for (size_type r= 0; r < A.dim1(); ++r) {
	MValue s(0);
	for (size_type k= r, i= 0; i < slots; ++i, k+= stride)
	    s+= A.ref_data()[k] * v[A.ref_minor()[k]];
	Assign::first_update(w[r], s);
    }
 }


// Row-major sparse_banded vector multiplication
template <typename MValue, typename MPara, typename VectorIn, typename VectorOut, typename Assign>
typename mtl::traits::enable_if_scalar<typename Collection<VectorOut>::value_type>::type
inline smat_cvec_mult(const sparse_banded<MValue, MPara>& A, const VectorIn& v, VectorOut& w, Assign, tag::row_major)
{
    vampir_trace<3069> tracer;
    typedef sparse_banded<MValue, MPara>                      Matrix;
    typedef typename Collection<VectorOut>::value_type        value_type;
    typedef typename Matrix::band_size_type                   band_size_type;
    typedef typename MPara::size_type                         size_type;
    typedef mtl::vec::dense_vector<band_size_type, vec::parameters<> > vector_type;

    if (size(w) == 0) return;
    const value_type z(math::zero(w[0]));

    size_type nr= num_rows(A), nc= num_cols(A), nb= A.ref_bands().size();
    if (nb == size_type(0) && Assign::init_to_zero) {
	set_to_zero(w);
	return;
    }

    vector_type bands(A.ref_bands()), begin_rows(max(0, -bands)), end_rows(min(nr, nc - bands));
    assert(end_rows[nb-1] > 0);

    // std::cout << "bands = " << bands << ", begin_rows = " << begin_rows << ", end_rows = " << end_rows << "\n";
    size_type begin_pos= 0, end_pos= nb - 1;

    // find lowest diagonal in row 0
    while (begin_pos < nb && begin_rows[begin_pos] > 0) begin_pos++;
    // if at the end, the first rows are empty
    if (begin_pos == nb && Assign::init_to_zero) {
	w[irange(begin_rows[--begin_pos])]= z;
	// std::cout << "w[0.." << begin_rows[begin_pos] << "] <- 0\n";
    }

    band_size_type from= begin_rows[begin_pos];
    // find first entry with same value
    while (begin_pos > 0 && begin_rows[begin_pos - 1] == from) {
	assert(from == 0); // should only happen when multiple bands start in row 0
	begin_pos--;
    }
    for (bool active= true; active; ) {
	// search backwards for the next-largest entry
	band_size_type to= begin_pos > 0 && begin_rows[begin_pos - 1] <= end_rows[end_pos] ? begin_rows[begin_pos - 1] : end_rows[end_pos];

	// std::cout << "rows " << from << ".." << to << ": with bands ";
	// for (size_type i= begin_pos; i <= end_pos; i++)
	//     std::cout << bands[i] << (i < end_pos ? ", " : "\n");

	const MValue* Aps= A.ref_data() + (from * nb + begin_pos);

	const band_size_type blocked_to= ((to - from) & -4) + from; 
	assert((blocked_to - from) % 4 == 0 && blocked_to >= band_size_type(from) && blocked_to <= band_size_type(to));
	for (band_size_type r= from; r < blocked_to; r+= 4) {
	    value_type     tmp0(z), tmp1(z), tmp2(z), tmp3(z);
	    const MValue   *Ap0= Aps, *Ap1= Aps + nb, *Ap2= Ap1 + nb, *Ap3= Ap2 + nb;
	    for (size_type b= begin_pos; b <= end_pos; ++b, ++Ap0, ++Ap1, ++Ap2, ++Ap3) {
		tmp0+= *Ap0 * v[r + bands[b]];
		tmp1+= *Ap1 * v[r + bands[b] + 1];
		tmp2+= *Ap2 * v[r + bands[b] + 2];
		tmp3+= *Ap3 * v[r + bands[b] + 3];
	    }
	    Assign::first_update(w[r], tmp0);
	    Assign::first_update(w[r+1], tmp1);
	    Assign::first_update(w[r+2], tmp2);
	    Assign::first_update(w[r+3], tmp3);
	    Aps+= 4 * nb;
	}

	for (band_size_type r= blocked_to; r < band_size_type(to); r++) {
	    value_type     tmp(z);
	    const MValue*  Ap= Aps;
	    for (size_type b= begin_pos; b <= end_pos; ++b, ++Ap)
		tmp+= *Ap * v[r + bands[b]];
	    Assign::first_update(w[r], tmp);
	    Aps+= nb;
	}
    
	if (begin_pos > 0) {
	    if (begin_rows[begin_pos-1] == to)
		begin_pos--;
	    if (end_rows[end_pos] == to)
		end_pos--;
	} else { // begin == 0 -> decrement end_pos or finish
	    if (end_rows[0] == to) {
		active= false;
		assert(end_rows[0] = end_rows[end_pos]);
	    } else {
		assert(end_pos > 0);
		end_pos--;
	    }
	}
	assert(begin_pos <= end_pos);
	from= to;
    }

    if (size_type(end_rows[0]) < nr  && Assign::init_to_zero) {
	w[irange(end_rows[0], nr)]= z;
	// std::cout << "w[" << end_rows[0] << ".." << nr << "] <- 0\n";
    }
}

// Sparse column-major matrix vector multiplication
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void smat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::col_major)
{
    vampir_trace<3023> tracer;
    using namespace tag; namespace traits = mtl::traits;
    using traits::range_generator;  
    using mtl::vec::set_to_zero;
    typedef typename range_generator<col, Matrix>::type       a_cur_type;             
    typedef typename range_generator<nz, a_cur_type>::type    a_icur_type;            

    typename traits::row<Matrix>::type                        row_a(A); 
    typename traits::const_value<Matrix>::type                value_a(A); 

    if (Assign::init_to_zero) set_to_zero(w);

    unsigned rv= 0; // traverse all rows of v
    for (a_cur_type ac= begin<col>(A), aend= end<col>(A); ac != aend; ++ac, ++rv) {
	typename Collection<VectorIn>::value_type    vv= v[rv]; 
	for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) 
	    Assign::update(w[row_a(*aic)], value_a(*aic) * vv);
    }
}

// Sparse matrix vector multiplication
template <typename Matrix, typename VectorIn, typename VectorOut, typename Assign>
inline void mat_cvec_mult(const Matrix& A, const VectorIn& v, VectorOut& w, Assign, tag::flat<tag::sparse>)
{
    smat_cvec_mult(A, v, w, Assign(), typename OrientedCollection<Matrix>::orientation());
}



}} // namespace mtl::matrix




#endif // MTL_MAT_VEC_MULT_INCLUDE

