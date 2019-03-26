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

#ifndef MTL_VECTOR_VEC_VEC_PMOP_EXPR_INCLUDE
#define MTL_VECTOR_VEC_VEC_PMOP_EXPR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/vec_expr.hpp>
#include <boost/numeric/mtl/operation/compute_summand.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>

namespace mtl { namespace vec {

// Model of VectorExpression
template <class E1, class E2, typename SFunctor>
struct vec_vec_pmop_expr
  : vec_expr< vec_vec_pmop_expr<E1, E2, SFunctor> >
{
    typedef vec_vec_pmop_expr<E1, E2, SFunctor>                   self;
    typedef typename mtl::operation::compute_summand<E1>::type    first_argument_type;
    typedef typename mtl::operation::compute_summand<E2>::type    second_argument_type;
    typedef SFunctor                                              functor_type;

    typedef typename SFunctor::result_type                        const_dereference_type;

    typedef const_dereference_type                                value_type;
    typedef typename Collection<first_argument_type>::size_type   size_type;

    vec_vec_pmop_expr( E1 const& v1, E2 const& v2 ) 
      : first(v1), second(v2)
    {
	first.value.delay_assign(); second.value.delay_assign();
    }

    void delay_assign() const {}

#if 0
    friend size_type inline size(const self& x)
    {
	assert( size(x.first.value) == 0 || size(x.first.value) == size(x.second.value) );
	return size(x.first.value);
    }
#endif

    template <typename EE1, typename EE2, typename SSFunctor>
    friend std::size_t size(const vec_vec_pmop_expr<EE1, EE2, SSFunctor>&);

    const_dereference_type operator() (size_type i) const
    {
        return SFunctor::apply(first.value(i), second.value(i));
    }

    const_dereference_type operator[] (size_type i) const
    {
        return SFunctor::apply(first.value(i), second.value(i));
    }

    operation::compute_summand<E1> const& first_argument() const { return first; }
    operation::compute_summand<E2> const& second_argument() const { return second; }

  private:
    operation::compute_summand<E1> first;
    operation::compute_summand<E2> second;
};

template <typename E1, typename E2, typename SFunctor>
inline std::size_t size(const vec_vec_pmop_expr<E1, E2, SFunctor>& v)
{
    MTL_CRASH_IF(mtl::vec::size(v.first.value) != 0 && mtl::vec::size(v.first.value) != mtl::vec::size(v.second.value),
	      "Incompatible size!");
    return mtl::vec::size(v.first.value);
}



}} // namespace mtl::vector

#endif // MTL_VECTOR_VEC_VEC_PMOP_EXPR_INCLUDE
