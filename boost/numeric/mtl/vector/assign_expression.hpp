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


#ifndef MTL_VEC_VEC_ASGN_EXPR_INCLUDE
#define MTL_VEC_VEC_ASGN_EXPR_INCLUDE


namespace mtl { namespace vec {

// Model of VectorExpression
template <class E1, class E2>
class vec_vec_asgn_expr 
{
public:
    // temporary solution
    typedef typename E1::value_type              value_type;
    // typedef typename glas::value< glas::scalar::vec_vec_asgn_expr<typename E1::value_type, typename E2::value_type > >::type value_type ;
    
    // temporary solution
    typedef typename E1::size_type               size_type;
    //typedef typename utilities::smallest< typename E1::size_type, typename E2::size_type >::type                          size_type ;

    typedef value_type const_dereference_type ;

    typedef E1 first_argument_type ;
    typedef E2 second_argument_type ;
    

    vec_vec_asgn_expr( first_argument_type& v1, second_argument_type const& v2 )
	: first( v1 ), second( v2 ), delayed_assign( false )
    {}

    ~vec_vec_asgn_expr()
    {
	if (!delayed_assign)
	    for (size_type i= 0; i < size(first); ++i)
		first( i )= second( i );
    }
    
    void delay_assign() { delayed_assign= true; }

    size_type size() const {
	assert( size(first) == second_.size() ) ;
	return size(first) ;
    }

     const_dereference_type operator() ( size_type i ) const {
	assert( delayed_assign );
        return first( i )= second( i ) ;
     }

     const_dereference_type operator[] ( size_type i ) const {
	assert( delayed_assign );
        return first( i )= second( i ) ;
     }

  private:
     first_argument_type&        first ;
     second_argument_type const& second ;
     bool                        delayed_assign;
  } ; // vec_vec_asgn_expr

} } // Namespace mtl::vector


namespace mtl { namespace traits {

  template <class E1, class E2>
  struct category< vec_vec_asgn_expr<E1,E2> > 
  {
      typedef vec type ;
  } ;

}} // Namespace mtl::traits

#endif

