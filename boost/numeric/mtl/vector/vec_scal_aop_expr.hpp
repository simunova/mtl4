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

#ifndef MTL_VEC_SCAL_AOP_EXPR_INCLUDE
#define MTL_VEC_SCAL_AOP_EXPR_INCLUDE

#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/vector/vec_expr.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace vec {

// Generic assign operation expression template for vectors
// Model of VectorExpression
template <class E1, class E2, typename SFunctor>
struct vec_scal_aop_expr 
    : public vec_expr< vec_scal_aop_expr<E1, E2, SFunctor> >
{
    typedef vec_expr< vec_scal_aop_expr<E1, E2, SFunctor> >  expr_base;
    typedef vec_scal_aop_expr                                self;
    typedef typename E1::value_type              value_type;
    
    // temporary solution
    typedef typename E1::size_type               size_type;
    //typedef typename utilities::smallest< typename E1::size_type, typename E2::size_type >::type                          size_type ;

    typedef value_type reference_type ;

    typedef E1 first_argument_type ;
    typedef E2 second_argument_type ;
    
    vec_scal_aop_expr( first_argument_type& v1, second_argument_type const& v2, bool delay= false )
      : first( v1 ), second( v2 ), delayed_assign( delay ), with_comma( false ), index(0)
    {}

    ~vec_scal_aop_expr()
    {
	if (!delayed_assign) {
	    vampir_trace<2018> tracer;
	    if (with_comma) {
		MTL_DEBUG_THROW_IF(index != mtl::vec::size(first), incompatible_size("Not all vector entries initialized!"));
	    } else
		for (size_type i= 0; i < mtl::vec::size(first); ++i)
		    SFunctor::apply( first(i), second );
	}
    }
    
    void delay_assign() const 
    { 
	MTL_DEBUG_THROW_IF(with_comma, logic_error("Comma notation conflicts with rich expression templates."));
	delayed_assign= true; 
    }

    //friend size_type inline size(const self& v) { return size(v.first); }
    template <typename EE1, typename EE2, typename SSFunctor>
    friend std::size_t size(const vec_scal_aop_expr<EE1, EE2, SSFunctor>& v);

    value_type& operator() ( size_type i ) const 
    {
	assert( delayed_assign && !with_comma);
	return SFunctor::apply( first(i), second );
    }

    value_type& operator[] ( size_type i ) const
    {
	assert( delayed_assign );
	return SFunctor::apply( first(i), second );
    }

    template <unsigned Offset>
    value_type& at(size_type i) const { 
	assert(delayed_assign);
	return SFunctor::apply(first(i+Offset), second);
    }

    template <typename Source>
    self& operator, (Source val)
    {
	//std::cout << "vec_scal_aop_expr::operator,\n";
	if (!with_comma) {
	    with_comma= true;
	    assert(index == 0);
	    SFunctor::apply( first(index++), second); // We haven't set v[0] yet
	}
	MTL_DEBUG_THROW_IF(index >= mtl::vec::size(first), range_error());
	SFunctor::apply( first(index++), val);
	return *this;
    }

  private:
    first_argument_type&         first ;
    second_argument_type const&  second ;
    mutable bool                 delayed_assign, with_comma;
    size_type                    index;
} ; // vec_scal_aop_expr


template <typename E1, typename E2, typename SFunctor>
inline std::size_t size(const vec_scal_aop_expr<E1, E2, SFunctor>& v)
{
    return size(v.first);
}



} } // Namespace mtl::vector





#endif

