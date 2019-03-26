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

#ifndef MTL_VECTOR_RVEC_MAT_MULT_INCLUDE
#define MTL_VECTOR_RVEC_MAT_MULT_INCLUDE

#include <cassert>
#include <boost/mpl/bool.hpp>

#include <boost/numeric/mtl/mtl_conditional_fwd.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>
#include <boost/numeric/mtl/operation/set_to_zero.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/operation/conj.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/meta_math/loop.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace vec {

namespace impl {

    template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1, typename Assign>
    struct fully_unroll_rvec_mat_mult
      : public meta_math::loop2<Index0, Max0, Index1, Max1>
    {
	typedef meta_math::loop2<Index0, Max0, Index1, Max1>                              base;
	typedef fully_unroll_rvec_mat_mult<base::next_index0, Max0, base::next_index1, Max1, Assign>  next_t;

	template <typename VectorIn, typename Matrix, typename VectorOut>
	static inline void apply(const VectorIn& v, const Matrix& A, VectorOut& w)
	{
	    Assign::update(w[base::index0], v[base::index1] * A[base::index1][base::index0]);
	    next_t::apply(v, A, w);
	}   
    };

    template <std::size_t Max0, std::size_t Max1, typename Assign>
    struct fully_unroll_rvec_mat_mult<Max0, Max0, Max1, Max1, Assign>
      : public meta_math::loop2<Max0, Max0, Max1, Max1>
    {
	typedef meta_math::loop2<Max0, Max0, Max1, Max1>                              base;

	template <typename VectorIn, typename Matrix, typename VectorOut>
	static inline void apply(const VectorIn& v, const Matrix& A, VectorOut& w)
	{
	    Assign::update(w[base::index0], v[base::index1] * A[base::index1][base::index0]);
	}   
    };

    struct rvec_mat_noop
    {
	template <typename VectorIn, typename Matrix, typename VectorOut>
	static inline void apply(const VectorIn&, const Matrix&, VectorOut&) {}
    };
} // impl

// Dense matrix vector multiplication with run-time matrix size
template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void dense_rvec_mat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, boost::mpl::true_)
{
    vampir_trace<1008> tracer;
    typedef typename static_num_rows<Matrix>::type size_type;
    static const size_type rows_a= static_num_rows<Matrix>::value, cols_a= static_num_cols<Matrix>::value;

    assert(rows_a > 0 && cols_a > 0);
    // w= A[all][0] * v[0];  N.B.: 1D is unrolled by the compiler faster (at least on gcc)
    for (size_type i= 0; i < cols_a; i++) 
	Assign::first_update(w[i], v[0] * A[0][i]);
	
    // corresponds to w+= v[1:] * A[1:][all]; if necessary
    typedef impl::fully_unroll_rvec_mat_mult<1, cols_a, 2, rows_a, Assign>  f2;
    typedef typename boost::mpl::if_c<(rows_a > 1), f2, impl::rvec_mat_noop>::type   f3;
    f3::apply(v, A, w);
}

// Dense matrix vector multiplication with run-time matrix size
template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void dense_rvec_mat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, boost::mpl::false_)
{
    // Naive implementation, will be moved to a functor and complemented with more efficient ones
	vampir_trace<3027> tracer;
    using math::zero; using mtl::vec::set_to_zero;
    if (size(w) == 0) return;

    if (Assign::init_to_zero) set_to_zero(w);

    typedef typename Collection<VectorOut>::value_type value_type;
    typedef typename Collection<Matrix>::size_type     size_type;

    for (size_type i= 0; i < num_cols(A); i++) {
	value_type tmp= zero(w[i]);
	for (size_type j= 0; j < num_rows(A); j++) 
	    tmp+= v[j] * A[j][i];
	Assign::update(w[i], tmp);
    }
}

// Dense vector matrix multiplication
template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void rvec_mat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, tag::flat<tag::dense>)
{
# ifdef MTL_NOT_UNROLL_FSIZE_MAT_VEC_MULT
    boost::mpl::false_        selector;
# else
    mtl::traits::is_static<Matrix> selector;
# endif
    dense_rvec_mat_mult(v, A, w, Assign(), selector);
}

// Multi-vector vector multiplication (tag::multi_vector is derived from dense)
template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void rvec_mat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, tag::flat<tag::multi_vector>)
{
	vampir_trace<2025> tracer;
    if (Assign::init_to_zero) set_to_zero(w);
    for (unsigned i= 0; i < num_cols(A); i++)
	Assign::update(w[i], dot_real(v, A.vector(i)));
}

// Transposed multi-vector vector multiplication (tag::transposed_multi_vector is derived from dense)
template <typename VectorIn, typename TransposedMatrix, typename VectorOut, typename Assign>
inline void rvec_mat_mult(const VectorIn& v, const TransposedMatrix& A, VectorOut& w, Assign, tag::flat<tag::transposed_multi_vector>)
{
	vampir_trace<2026> tracer;
    typename TransposedMatrix::const_ref_type B= A.ref; // Referred matrix

    if (Assign::init_to_zero) set_to_zero(w);
    for (unsigned i= 0; i < num_cols(B); i++)
	Assign::update(w, v[i] * trans(B.vector(i)));
}

// Hermitian multi-vector vector multiplication (tag::hermitian_multi_vector is derived from dense)
template <typename VectorIn, typename HermitianMatrix, typename VectorOut, typename Assign>
inline void rvec_mat_mult(const VectorIn& v, const HermitianMatrix& A, VectorOut& w, Assign, tag::flat<tag::hermitian_multi_vector>)
{
	vampir_trace<2027> tracer;
    typename HermitianMatrix::const_ref_type B= A.const_ref(); // Referred matrix

    if (Assign::init_to_zero) set_to_zero(w);
    for (unsigned i= 0; i < num_cols(B); i++)
	Assign::update(w, v[i] * mtl::vec::conj(trans(B.vector(i))));
}


// Vector sparse matrix multiplication
template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void rvec_mat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, tag::flat<tag::sparse>)
{
	vampir_trace<3028> tracer;
    rvec_smat_mult(v, A, w, Assign(), typename OrientedCollection<Matrix>::orientation());
}



template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void rvec_smat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, tag::row_major)
{
	using namespace tag; namespace traits = mtl::traits;
	using traits::range_generator;  
	using mtl::vec::set_to_zero;
        typedef typename range_generator<row, Matrix>::type       a_cur_type;             
        typedef typename range_generator<nz, a_cur_type>::type    a_icur_type;            

	typename traits::col<Matrix>::type                        col_a(A); 
	typename traits::const_value<Matrix>::type                value_a(A); 

	if (Assign::init_to_zero) set_to_zero(w);

	unsigned cv= 0; // traverse all columns of v
	for (a_cur_type ac= begin<row>(A), aend= end<row>(A); ac != aend; ++ac, ++cv) {
	    typename Collection<VectorIn>::value_type    vv= v[cv]; 
	    for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) 
		Assign::update(w[col_a(*aic)], vv * value_a(*aic));
	}
}

template <typename VectorIn, typename Matrix, typename VectorOut, typename Assign>
inline void rvec_smat_mult(const VectorIn& v, const Matrix& A, VectorOut& w, Assign, tag::col_major)
{
    using namespace tag; 
    using mtl::traits::range_generator;  
    using math::zero;
    using mtl::vec::set_to_zero;

    typedef typename range_generator<col, Matrix>::type       a_cur_type;    
    typedef typename range_generator<nz, a_cur_type>::type    a_icur_type;            
    typename mtl::traits::row<Matrix>::type                   row_a(A); 
    typename mtl::traits::const_value<Matrix>::type           value_a(A); 

    if (Assign::init_to_zero) set_to_zero(w);

    typedef typename Collection<VectorOut>::value_type        value_type;

    a_cur_type ac= begin<col>(A), aend= end<col>(A);
    for (unsigned i= 0; ac != aend; ++ac, ++i) {
	value_type tmp= zero(w[i]);
	for (a_icur_type aic= begin<nz>(ac), aiend= end<nz>(ac); aic != aiend; ++aic) 
	    tmp+= v[row_a(*aic)] * value_a(*aic);	
	Assign::update(w[i], tmp);
    }
}

}} // namespace mtl::vector

#endif // MTL_VECTOR_RVEC_MAT_MULT_INCLUDE
