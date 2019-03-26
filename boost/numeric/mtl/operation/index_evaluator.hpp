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

#ifndef MTL_INDEX_EVALUATOR_INCLUDE
#define MTL_INDEX_EVALUATOR_INCLUDE

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/operation/lazy_assign.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/vector/vec_vec_aop_expr.hpp>
#include <boost/numeric/mtl/vector/vec_scal_aop_expr.hpp>
#include <boost/numeric/mtl/vector/reduction_index_evaluator.hpp>
#include <boost/numeric/mtl/vector/dot_index_evaluator.hpp>
#include <boost/numeric/mtl/vector/row_mat_cvec_index_evaluator.hpp>
#include <boost/numeric/mtl/operation/fused_index_evaluator.hpp>
#include <boost/numeric/mtl/operation/fused_expr.hpp>
#include <boost/numeric/itl/pc/ic_0.hpp>
#include <boost/numeric/itl/pc/ilu_0.hpp>

namespace mtl {

/// Overloaded function that maps from lazy expressions to the according index-wise evaluation classes
template <typename T, typename U, typename Assign>
typename boost::enable_if<boost::mpl::and_<mtl::traits::is_vector<T>, mtl::traits::is_vector<U> >, 
			  mtl::vec::vec_vec_aop_expr<T, U, Assign> >::type
inline index_evaluator(lazy_assign<T, U, Assign>& lazy)
{
    return mtl::vec::vec_vec_aop_expr<T, U, Assign>(lazy.first, lazy.second, true);
}

template <typename T, typename U, typename Assign>
typename boost::enable_if<boost::mpl::and_<mtl::traits::is_vector<T>, mtl::traits::is_scalar<U> >, 
			  mtl::vec::vec_scal_aop_expr<T, U, Assign> >::type
inline index_evaluator(lazy_assign<T, U, Assign>& lazy)
{
    return mtl::vec::vec_scal_aop_expr<T, U, Assign>(lazy.first, lazy.second, true);
}

template <typename Scalar, typename Vector, typename Functor, typename Assign>
mtl::vec::reduction_index_evaluator<Scalar, Vector, Functor, Assign>
inline index_evaluator(lazy_assign<Scalar, mtl::vec::lazy_reduction<Vector, Functor>, Assign>& lazy)
{
    return mtl::vec::reduction_index_evaluator<Scalar, Vector, Functor, Assign>(lazy.first, lazy.second.v);
}

template <typename Scalar, unsigned long Unroll, typename Vector1, 
	  typename Vector2, typename ConjOpt, typename Assign>
mtl::vec::dot_index_evaluator<Scalar, Vector1, Vector2, ConjOpt, Assign>
inline index_evaluator(lazy_assign<Scalar, mtl::vec::dot_class<Unroll, Vector1, Vector2, ConjOpt>, Assign>& lazy)
{
    return mtl::vec::dot_index_evaluator<Scalar, Vector1, Vector2, ConjOpt, Assign>(lazy.first, lazy.second.v1, lazy.second.v2);
}

template <typename VectorOut, typename Matrix, typename VectorIn, typename Assign>
mtl::vec::row_mat_cvec_index_evaluator<VectorOut, Matrix, VectorIn, Assign>
inline index_evaluator(lazy_assign<VectorOut, mtl::mat_cvec_times_expr<Matrix, VectorIn>, Assign>& lazy)
{
    return mtl::vec::row_mat_cvec_index_evaluator<VectorOut, Matrix, VectorIn, Assign>(lazy.first, lazy.second.first, lazy.second.second);
}

template <typename V1, typename Matrix, typename Value, typename V2>
itl::pc::ic_0_evaluator<V1, itl::pc::solver<itl::pc::ic_0<Matrix, Value>, V2, true> >
inline index_evaluator(lazy_assign<V1, itl::pc::solver<itl::pc::ic_0<Matrix, Value>, V2, true>, assign::assign_sum>& lazy)
{
    return itl::pc::ic_0_evaluator<V1, itl::pc::solver<itl::pc::ic_0<Matrix, Value>, V2, true> >(lazy.first, lazy.second);
}

// reuse ic_0_evaluator with ilu_0 preconditioner since IC(0) and ILU(0) do the same at the upper triangle ;-)
template <typename V1, typename Factorizer, typename Matrix, typename Value, typename V2>
itl::pc::ic_0_evaluator<V1, itl::pc::solver<itl::pc::ilu<Matrix, Factorizer, Value>, V2, true> >
inline index_evaluator(lazy_assign<V1, itl::pc::solver<itl::pc::ilu<Matrix, Factorizer, Value>, V2, true>, assign::assign_sum>& lazy)
{
    return itl::pc::ic_0_evaluator<V1, itl::pc::solver<itl::pc::ilu<Matrix, Factorizer, Value>, V2, true> >(lazy.first, lazy.second);
}

template <typename T, typename U>
inline mtl::vec::fused_index_evaluator<T, U> index_evaluator(fused_expr<T, U>& expr)
{  return mtl::vec::fused_index_evaluator<T, U>(expr.first, expr.second); }

} // namespace mtl

#endif // MTL_INDEX_EVALUATOR_INCLUDE
