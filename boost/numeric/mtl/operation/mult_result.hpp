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

#ifndef MTL_MULT_RESULT_INCLUDE
#define MTL_MULT_RESULT_INCLUDE

#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

#if 0
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/matrix/mat_mat_times_expr.hpp>
#include <boost/numeric/mtl/matrix/implicit_dense.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/operation/mat_cvec_times_expr.hpp>
#include <boost/numeric/mtl/vector/rvec_mat_times_expr.hpp>
#endif

namespace mtl { namespace traits {

template <typename Op1, typename Op2, typename MultOp> struct mult_result_aux;
template <typename Op1, typename Op2, typename MultOp> struct vec_mult_result_aux;
    //template <typename Op1, typename Op2, typename MultOp1, typename MultOp2> struct mult_result_if_equal_aux;

/// Result type for multiplying arguments of types Op1 and Op2
/** Can be used in enable-if-style as type is only defined when appropriate. 
    This one is used if at least one argument is a matrix.
**/
template <typename Op1, typename Op2>
struct mult_result 
  : public mult_result_aux<Op1, Op2, typename ashape::mult_op<typename ashape::ashape<Op1>::type, 
							      typename ashape::ashape<Op2>::type >::type>
{}; 


/// Result type for multiplying arguments of types Op1 and Op2
/** Can be used in enable-if-style as type is only defined when appropriate. 
    This one is used if at least one argument is a vector and none is a matrix.
**/
template <typename Op1, typename Op2>
struct vec_mult_result 
  : public vec_mult_result_aux<Op1, Op2, typename ashape::mult_op<typename ashape::ashape<Op1>::type, 
								  typename ashape::ashape<Op2>::type >::type>
{}; 


/// Result type for multiplying arguments of types Op1 and Op2
/** MultOp according to the algebraic shapes **/
template <typename Op1, typename Op2, typename MultOp>
struct mult_result_aux {};

/// Scale matrix from left
template <typename Op1, typename Op2>
struct mult_result_aux<Op1, Op2, ::mtl::ashape::scal_mat_mult> 
{
    typedef mat::scaled_view<Op1, Op2> type;
};


/// Scale matrix from right needs functor for scaling from right
template <typename Op1, typename Op2>
struct mult_result_aux<Op1, Op2, ::mtl::ashape::mat_scal_mult> 
{
    typedef mat::rscaled_view<Op1, Op2> type;
};

/// Multiply matrices
template <typename Op1, typename Op2>
struct mult_result_aux<Op1, Op2, ::mtl::ashape::mat_mat_mult> 
{
    typedef mat::mat_mat_times_expr<Op1, Op2> type;
};

/// Multiply matrix with column vector
template <typename Op1, typename Op2>
struct mult_result_aux<Op1, Op2, ::mtl::ashape::mat_cvec_mult> 
{
    typedef mat_cvec_times_expr<Op1, Op2> type;
};

/// Multiply column with row vector and return implicit matrix
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::cvec_rvec_mult> 
{
    typedef mtl::mat::outer_product_matrix<Op1, Op2> type;
};

/// Result type for multiplying arguments of types Op1 and Op2
/** MultOp according to the algebraic shapes **/
template <typename Op1, typename Op2, typename MultOp>
struct vec_mult_result_aux {};

/// Scale row vector from left
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::scal_rvec_mult> 
{
    typedef vec::scaled_view<Op1, Op2> type;
};




/// Scale column vector from left
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::scal_cvec_mult> 
{
    typedef vec::scaled_view<Op1, Op2> type;
};

/// Multiply row vector with matrix
// added by Cornelius and Peter
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::rvec_mat_mult> 
{
    typedef vec::rvec_mat_times_expr<Op1, Op2> type;
};

/// Scale row vector from right
// added by Hui Li
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::rvec_scal_mult> 
{
    typedef vec::rscaled_view<Op1, Op2> type;
};

/// Scale column vector from right
// added by Hui Li
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::cvec_scal_mult> 
{
    typedef vec::rscaled_view<Op1, Op2> type;
};
	

/// Enabler if operation is rvec_cvec_mult
template <typename Op1, typename Op2, typename Result>
struct lazy_enable_if_rvec_cvec_mult 
  : boost::lazy_enable_if<boost::is_same<typename ashape::mult_op<typename ashape::ashape<Op1>::type, 
								  typename ashape::ashape<Op2>::type >::type,
					 ashape::rvec_cvec_mult>,
			  Result>
{};

#ifndef MTL_WITHOUT_VECTOR_ELE_OPS
/// Product of row vectors is performed element-wise (unless disabled with MTL_WITHOUT_VECTOR_ELE_OPS)
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::rvec_rvec_mult> 
{
    typedef typename Collection<Op1>::value_type v1;
    typedef typename Collection<Op2>::value_type v2;

    typedef vec::vec_vec_op_expr<Op1, Op2, mtl::sfunctor::times<v1, v2> > type;
};

/// Product of column vectors is performed element-wise (unless disabled with MTL_WITHOUT_VECTOR_ELE_OPS)
template <typename Op1, typename Op2>
struct vec_mult_result_aux<Op1, Op2, ::mtl::ashape::cvec_cvec_mult> 
  : vec_mult_result_aux<Op1, Op2, ::mtl::ashape::rvec_rvec_mult> 
{};
#endif // MTL_WITHOUT_VECTOR_ELE_OPS


}} // namespace mtl::traits

#endif // MTL_MULT_RESULT_INCLUDE
