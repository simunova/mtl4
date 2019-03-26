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

#ifndef MTL_ROW_MAT_CVEC_INDEX_EVALUATOR_INCLUDE
#define MTL_ROW_MAT_CVEC_INDEX_EVALUATOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>

namespace mtl { namespace vec {

/// ET class to evaluate the product of a row-major matrix and a column vector row by row
template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign>
struct row_mat_cvec_index_evaluator
{
    MTL_STATIC_ASSERT((mtl::traits::is_row_major<Matrix>::value), "Only row-major matrices supported here.");
    typedef typename mtl::Collection<VectorOut>::value_type        value_type;
    typedef typename mtl::Collection<Matrix>::size_type            size_type; 

    row_mat_cvec_index_evaluator(VectorOut& w, const Matrix& A, const VectorIn& v) : w(w), A(A), v(v) {}

    template <unsigned Offset>
    void at(size_type i, boost::mpl::true_)
    {
	// value_type tmp(math::zero(w[i+Offset]));
	value_type tmp(0);
	const size_type cj0= A.ref_major()[i+Offset], cj1= A.ref_major()[i+Offset+1];
	for (size_type j= cj0; j != cj1; ++j)
	    tmp+= A.data[j] * v[A.ref_minor()[j]];
	Assign::first_update(w[i+Offset], tmp);
    }

    template <unsigned Offset>
    void at(size_type i, boost::mpl::false_)
    {
	value_type tmp(math::zero(w[i+Offset]));
	for (size_type j= 0; j < num_cols(A); j++) 
	    tmp+= A[i][j] * v[j];
	Assign::first_update(w[i+Offset], tmp);
    }

    template <unsigned Offset>
    void at(size_type i)
    { 
	at<Offset>(i, typename mtl::traits::is_sparse<Matrix>::type());
    }

    void operator()(size_type i) { at<0>(i); }
    void operator[](size_type i) { at<0>(i); }

    VectorOut&      w;
    const Matrix&   A;
    const VectorIn& v;
};

template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign>
inline std::size_t size(const row_mat_cvec_index_evaluator<VectorOut, Matrix, VectorIn, Assign>& eval)
{
    return size(eval.w);
}


}} // namespace mtl::vector

#endif // MTL_ROW_MAT_CVEC_INDEX_EVALUATOR_INCLUDE
