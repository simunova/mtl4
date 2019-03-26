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

#ifndef MTL_MAT_MAT_ELE_TIMES_EXPR_INCLUDE
#define MTL_MAT_MAT_ELE_TIMES_EXPR_INCLUDE

#include <boost/shared_ptr.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/mat_mat_op_expr.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>


namespace mtl { namespace mat {

template <typename E1, typename E2>
struct mat_mat_ele_times_expr 
    : public mat_mat_op_expr< E1, E2, mtl::sfunctor::times<typename E1::value_type, typename E2::value_type> >,
      public mat_expr< mat_mat_ele_times_expr<E1, E2> >
{
    typedef mat_mat_op_expr< E1, E2, mtl::sfunctor::times<typename E1::value_type, typename E2::value_type> > op_base;
    typedef mat_expr< mat_mat_ele_times_expr<E1, E2> >                                                        crtp_base;
    typedef mat_mat_ele_times_expr               self;
    typedef E1                                   first_argument_type ;
    typedef E2                                   second_argument_type ;
    typedef typename E1::orientation             orientation;
    typedef mtl::non_fixed::dimensions           dim_type;
    typedef typename E1::key_type                key_type;

    
    mat_mat_ele_times_expr( E1 const& v1, E2 const& v2 )
	: op_base( v1, v2 ), crtp_base(*this), first(v1), second(v2)
    {}

    first_argument_type const&  first ;
    second_argument_type const& second ;
};

template <typename E1, typename E2>
std::size_t inline num_rows(const mat_mat_ele_times_expr<E1, E2>& expr) 
{ return num_rows(expr.first); }

template <typename E1, typename E2>
std::size_t inline num_cols(const mat_mat_ele_times_expr<E1, E2>& expr) 
{ return num_cols(expr.second); }

template <typename E1, typename E2>
std::size_t inline size(const mat_mat_ele_times_expr<E1, E2>& expr) 
{ return num_rows(expr) * num_cols(expr); }

}} // Namespace mtl::matrix

#endif // MTL_MAT_MAT_ELE_TIMES_EXPR_INCLUDE
