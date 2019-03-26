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

// Adapted from GLAS implementation by Karl Meerbergen and Toon Knappen


#ifndef MTL_VEC_VEC_OP_EXPR_INCLUDE
#define MTL_VEC_VEC_OP_EXPR_INCLUDE

#include <boost/numeric/mtl/vector/vec_expr.hpp>

namespace mtl { namespace vec {

// Model of VectorExpression
template <typename E1, typename E2, typename SFunctor>
struct vec_vec_op_expr 
  : public vec_expr< vec_vec_op_expr<E1, E2, SFunctor> >
{
    typedef vec_expr< vec_vec_op_expr<E1, E2, SFunctor> > expr_base;
    typedef vec_vec_op_expr                     self;

    // temporary solution
    // typedef typename E1::value_type              value_type;
    // typedef typename glas::value< glas::scalar::vec_vec_op_expr<typename E1::value_type, typename E2::value_type > >::type value_type ;
    typedef typename SFunctor::result_type       value_type;

    // temporary solution
    typedef typename E1::size_type               size_type;
    //typedef typename utilities::smallest< typename E1::size_type, typename E2::size_type >::type                          size_type ;

    typedef value_type                           const_dereference_type ;

    typedef E1                                   first_argument_type ;
    typedef E2                                   second_argument_type ;
    
public:
    vec_vec_op_expr( first_argument_type const& v1, second_argument_type const& v2 )
	: first( v1 ), second( v2 )
    {
	// std::cerr << "vec_vec_op_expr.vec_vec_op_expr()  " << size(first) << "  " << size(second) << "\n";
	first.delay_assign();
	second.delay_assign();
    }
    
    void delay_assign() const {}

    template <typename EE1, typename EE2, typename SSFunctor>
    friend std::size_t size(const vec_vec_op_expr<EE1, EE2, SSFunctor>& v);

    const_dereference_type operator() ( size_type i ) const
    {
        return SFunctor::apply( first( i ), second( i ) ) ;
    }

    const_dereference_type operator[] ( size_type i ) const
    {
        return SFunctor::apply( first( i ), second( i ) ) ;
    }

  private:
    first_argument_type const&  first ;
    second_argument_type const& second ;
} ; // vec_vec_op_expr


template <typename E1, typename E2, typename SFunctor>
inline std::size_t size(const vec_vec_op_expr<E1, E2, SFunctor>& v)
{
    // std::cerr << "vec_vec_op_expr.size(): " << size(v.first) << " vs. " << size(v.second) << "\n";
    assert( size(v.first) == size(v.second) ) ;
    return size(v.first) ;
}

    
} } // Namespace mtl::vector


#endif

