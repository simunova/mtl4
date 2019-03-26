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

#ifndef MTL_RVEC_MAT_TIMES_EXPR_INCLUDE
#define MTL_RVEC_MAT_TIMES_EXPR_INCLUDE

#include <boost/numeric/mtl/operation/bin_op_expr.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace vec {

template <typename E1, typename E2>
struct rvec_mat_times_expr 
  : public bin_op_expr< E1, E2 >,
    public mtl::vec::vec_expr< rvec_mat_times_expr<E1, E2> >
{
    typedef bin_op_expr< E1, E2 >         base;
    typedef rvec_mat_times_expr<E1, E2>   self;

    typedef typename Multiplicable<typename Collection<E1>::value_type,
				   typename Collection<E2>::value_type>::result_type   value_type;
    typedef typename Collection<E2>::size_type    size_type;
  
    rvec_mat_times_expr( E1 const& v, E2 const& A ) : base(v, A) {}
};

template <typename E1, typename E2>
std::size_t inline size(const rvec_mat_times_expr<E1, E2>& x)
{ return num_cols(x.first); }

}} // namespace mtl::vector

#endif // MTL_RVEC_MAT_TIMES_EXPR_INCLUDE
