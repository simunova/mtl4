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

#ifndef MTL_EVALUATE_LAZY_INCLUDE
#define MTL_EVALUATE_LAZY_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/operation/lazy_assign.hpp>
#include <boost/numeric/mtl/vector/lazy_reduction.hpp>


namespace mtl {

/// Overloaded function to evaluate lazy expressions
template <typename T, typename U, typename Assign>
void inline evaluate_lazy(lazy_assign<T, U, Assign>& lazy)
{
    Assign::first_update(lazy.first, lazy.second);
}

template <typename T, typename Vector, typename Functor, typename Assign>
void inline evaluate_lazy(lazy_assign<T, vec::lazy_reduction<Vector, Functor>, Assign>& lazy)
{
    lazy.first= Functor::post_reduction(vec::reduction<4, Functor, T>::apply(lazy.second.v));
}

// To be checked whether this works
template <typename T, typename U> 
void inline evaluate_lazy(fused_expr<T,  U>& fu)
{
    evaluate_lazy(fu.first);
    evaluate_lazy(fu.second);
}



} // namespace mtl

#endif // MTL_EVALUATE_LAZY_INCLUDE
