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

#ifndef MTL__OPERATION_COMPUTE_SUMMAND_INCLUDE
#define MTL__OPERATION_COMPUTE_SUMMAND_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/copy_expression_const_ref_container.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace operation {

/// Compute a summand in an expression
/** For instance matrix vector products are transformed into mult function calls
    when assigned to a column vector. Adding such a vector to other vector expressions
    requires to compute (evaluate) the summand first and then add the resulting 
    vector. The current implementation assumes that the result of the operation
    can be represented by the multiplied column vector. This is for instance wrong
    when the matrix is complex and the vector real. To handle this is signifantly
    more complicated and is planned for the future.
**/
template <typename Expr>
struct compute_summand
{
    typedef Expr type;
    compute_summand(const Expr& expr) : value(expr) {}
    // value is a const& if Expr is a true vector and a copy if it is an expression
    typename mtl::traits::copy_expression_const_ref_container<Expr>::type value;
};


/// Specialization for matrix vector products
template <typename Matrix, typename CVector>
struct compute_summand< mat_cvec_times_expr<Matrix, CVector> >
{
    typedef CVector    type;

    compute_summand(const mat_cvec_times_expr<Matrix, CVector>& expr) 
      : value(num_rows(expr.first))
#ifndef NDEBUG // might be helpful, e.g. for generating AST
      , first(expr.first), second(expr.second)
#endif
    {
	vampir_trace<3005> tracer;
	value= expr.first * expr.second;
    }

    CVector value;
#ifndef NDEBUG
    const Matrix&  first;
    const CVector& second;
#endif
};
	
/// Specialization for matrix vector products
template <typename Matrix, typename CVector>
struct compute_summand< vec::mat_cvec_multiplier<Matrix, CVector> >
{
    typedef CVector    type;

    compute_summand(const vec::mat_cvec_multiplier<Matrix, CVector>& expr) 
      : value(num_rows(expr.A))
    {
	vampir_trace<3005> tracer;
	expr.assign_to(value);
    }

    CVector value;
};
	
}} // namespace mtl::operation

#endif // MTL__OPERATION_COMPUTE_SUMMAND_INCLUDE
