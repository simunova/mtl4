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

#ifndef MTL_MAT_MAT_OP_EXPR_INCLUDE
#define MTL_MAT_MAT_OP_EXPR_INCLUDE

#include <boost/numeric/mtl/matrix/mat_expr.hpp>
#include <boost/numeric/mtl/matrix/crtp_base_matrix.hpp>

namespace mtl { namespace mat {

    
template <typename E1, typename E2, typename SFunctor>
struct mat_mat_op_expr
    : public const_crtp_matrix_bracket< mat_mat_op_expr<E1, E2, SFunctor>, 
					typename SFunctor::result_type, 
					typename E1::size_type >
   // : public mat_expr< mat_mat_op_expr<E1, E2, SFunctor> >
{
    // typedef mat_expr< mat_mat_op_expr<E1, E2, SFunctor> > expr_base;
    typedef mat_mat_op_expr                               self;

    typedef typename SFunctor::result_type      value_type;

    // temporary solution
    typedef typename E1::size_type               size_type;

    typedef value_type                           const_dereference_type ;

    typedef E1                                   first_argument_type ;
    typedef E2                                   second_argument_type ;
    
    mat_mat_op_expr( first_argument_type const& v1, second_argument_type const& v2 )
	: // expr_base( *this ), 
	  first( v1 ), second( v2 )
    {
#if 0
	first.delay_assign();
	second.delay_assign();
#endif
    }
    
    void delay_assign() const {}

    void check_shape() const {} // consistency of shapes depend on operation

    const_dereference_type operator() (size_type row, size_type col) const
    {
        return SFunctor::apply( first(row, col), second(row, col) ) ;
	// return SFunctor::apply( first[row][col], second[row][col] ) ;
    }

#if 0
    template <typename> friend size_type size(const self&);
    template <typename> friend size_type num_rows(const self&);
    template <typename> friend size_type num_cols(const self&);
#endif

    first_argument_type const&  first ;
    second_argument_type const& second ;
};


template <typename E1, typename E2, typename SFunctor>
typename mat_mat_op_expr<E1, E2, SFunctor>::size_type 
inline size(mat_mat_op_expr<E1, E2, SFunctor> const& expr)
{
    expr.check_shape();
    return size(expr.first) ;
}

template <typename E1, typename E2, typename SFunctor>
typename mat_mat_op_expr<E1, E2, SFunctor>::size_type 
inline num_rows(mat_mat_op_expr<E1, E2, SFunctor> const& expr)
{
    expr.check_shape();
    return num_rows(expr.first) ;
}

template <typename E1, typename E2, typename SFunctor>
typename mat_mat_op_expr<E1, E2, SFunctor>::size_type 
inline num_cols(mat_mat_op_expr<E1, E2, SFunctor> const& expr)
{
    expr.check_shape();
    return num_cols(expr.first) ;
}


}} // namespace mtl

#endif // MTL_MAT_MAT_OP_EXPR_INCLUDE
