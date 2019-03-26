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

#ifndef MTL_MAT_MAT_PLUS_EXPR_INCLUDE
#define MTL_MAT_MAT_PLUS_EXPR_INCLUDE

#include <boost/numeric/mtl/matrix/mat_mat_op_expr.hpp>
#include <boost/numeric/mtl/vector/vec_vec_pmop_expr.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>

namespace mtl { namespace mat {

template <typename E1, typename E2>
struct mat_mat_plus_expr 
  : public mat_mat_op_expr< E1, E2, mtl::sfunctor::plus<typename E1::value_type, typename E2::value_type> >,
    public mat_expr< mat_mat_plus_expr<E1, E2> >
{
    typedef mat_mat_op_expr< E1, E2, mtl::sfunctor::plus<typename E1::value_type, typename E2::value_type> > op_base;
    typedef mat_expr< mat_mat_plus_expr<E1, E2> >                                                       crtp_base;
    typedef E1                                   first_argument_type ;
    typedef E2                                   second_argument_type ;
    
    mat_mat_plus_expr( E1 const& v1, E2 const& v2 )
      : op_base( v1, v2 ), crtp_base(*this), first(v1), second(v2)
    {}

    first_argument_type const&  first ;
    second_argument_type const& second ;
};

// Same as mat_mat_plus_expr for pair of dense matrix expressions
// Future versions will probably provide more efficient implementations for it
template <typename E1, typename E2>
struct dmat_dmat_plus_expr 
  : public mat_mat_op_expr< E1, E2, mtl::sfunctor::plus<typename E1::value_type, typename E2::value_type> >
{
    typedef mat_mat_op_expr< E1, E2, mtl::sfunctor::plus<typename E1::value_type, typename E2::value_type> > base;
    dmat_dmat_plus_expr( E1 const& v1, E2 const& v2 )
      : base( v1, v2 )
    {}
};
    
template <typename E1, typename E2>
struct mv_mv_plus_expr
  : mat_mat_plus_expr<E1, E2>
{
    typedef mat_mat_plus_expr< E1, E2 > base;
    typedef typename E1::vector_type V1;
    typedef typename E2::vector_type V2;
    typedef mtl::vec::vec_vec_pmop_expr< V1, V2, mtl::sfunctor::plus<typename V1::value_type, typename V2::value_type> > vector_type;

    mv_mv_plus_expr( E1 const& v1, E2 const& v2 )
      : base( v1, v2 )
    {}

    vector_type vector(std::size_t c) const { return vector_type(this->first.vector(c), this->second.vector(c)); }
};




}} // Namespace mtl::matrix

#endif // MTL_MAT_MAT_PLUS_EXPR_INCLUDE
