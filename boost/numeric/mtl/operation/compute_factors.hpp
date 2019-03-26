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

#ifndef MTL_COMPUTE_FACTORS_INCLUDE
#define MTL_COMPUTE_FACTORS_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace operation {

// Default is to just refer to the expression
template <typename Result, typename Expr>
struct compute_one_factor 
{
    typedef const Expr&           type;
    typedef const Expr&           const_reference;

    compute_one_factor(type src) : value(src) {}

    type value;
};

template <typename Result, typename E1, typename E2>
struct compute_one_factor<Result, mat::mat_mat_times_expr<E1, E2> > 
{
    typedef Result                type;
    typedef const Result&         const_reference;

    compute_one_factor(const mat::mat_mat_times_expr<E1, E2>& src) 
	: value(src.first * src.second) {}

    type value;
};

template <typename Result, typename E1, typename E2>
struct compute_one_factor<Result, mat::mat_mat_ele_times_expr<E1, E2> > 
{
    typedef Result                type;
    typedef const Result&         const_reference;

    compute_one_factor(const mat::mat_mat_ele_times_expr<E1, E2>& src) 
	: value(ele_prod(src.first, src.second)) {}

    type value;
};


// Only defined for mat::mat_mat_times_expr and mat::mat_mat_ele_times_expr
template <typename Result, typename Expr>
struct compute_factors {};

template <typename Result, typename E1, typename E2>
struct compute_factors<Result, mat::mat_mat_times_expr<E1, E2> >
{
    compute_factors(const mat::mat_mat_times_expr<E1, E2>& src) 
	: first_factor(src.first), second_factor(src.second),
	  first(first_factor.value), second(second_factor.value)
    {}

    compute_one_factor<Result, E1> first_factor;
    compute_one_factor<Result, E2> second_factor;

    typename compute_one_factor<Result, E1>::const_reference first;
    typename compute_one_factor<Result, E2>::const_reference second;
};

// First factor is handled implicitly in the evaluation
template <typename Result, typename E1, typename E2>
struct compute_factors<Result, mat::mat_mat_ele_times_expr<E1, E2> >
{
    compute_factors(const mat::mat_mat_ele_times_expr<E1, E2>& src) 
	: first(src.first), second_factor(src.second),
	  second(second_factor.value)
    {}

    const E1&                                                first;
    compute_one_factor<Result, E2>                           second_factor;
    typename compute_one_factor<Result, E2>::const_reference second;
};

}} // namespace mtl::operation

#endif // MTL_COMPUTE_FACTORS_INCLUDE
