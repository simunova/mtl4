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

#ifndef MTL_TRAITS_COPY_EXPRESSION_CONST_REF_CONTAINER_INCLUDE
#define MTL_TRAITS_COPY_EXPRESSION_CONST_REF_CONTAINER_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Type trait that defines the type itself for expressions and a const reference otherwise
/** Needed for re-constructed expression templates because references to expression are 
    often invalid. **/
template <typename Container>
struct copy_expression_const_ref_container
{
    typedef const Container& type;
};


template <class E1, class E2, typename SFunctor>
struct copy_expression_const_ref_container< mtl::vec::vec_vec_pmop_expr<E1, E2, SFunctor> >
{
    typedef mtl::vec::vec_vec_pmop_expr<E1, E2, SFunctor> type;
};

template <class E1, class E2, typename SFunctor>
struct copy_expression_const_ref_container< mtl::vec::vec_vec_aop_expr<E1, E2, SFunctor> >
{
    typedef mtl::vec::vec_vec_aop_expr<E1, E2, SFunctor> type;
};

template <typename Functor, typename Vector> 
struct copy_expression_const_ref_container< mtl::vec::map_view<Functor, Vector> >
{
    typedef mtl::vec::map_view<Functor, Vector> type;
};

template <typename Scaling, typename Vector>
struct copy_expression_const_ref_container< mtl::vec::scaled_view<Scaling, Vector> >
{
    typedef mtl::vec::scaled_view<Scaling, Vector> type;
};

template <typename Vector, typename RScaling>
struct copy_expression_const_ref_container< mtl::vec::rscaled_view<Vector, RScaling> >
{
    typedef mtl::vec::rscaled_view<Vector, RScaling> type;
};

template <typename Vector, typename Divisor>
struct copy_expression_const_ref_container< mtl::vec::divide_by_view<Vector, Divisor> >
{
    typedef mtl::vec::divide_by_view<Vector, Divisor> type;
};

template <typename Vector>
struct copy_expression_const_ref_container< mtl::vec::conj_view<Vector> >
{
    typedef mtl::vec::conj_view<Vector> type;
};





}} // namespace mtl::traits

#endif // MTL_TRAITS_COPY_EXPRESSION_CONST_REF_CONTAINER_INCLUDE
