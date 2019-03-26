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

#ifndef MTL_VEC_VEC_AOP_EXPR_INCLUDE
#define MTL_VEC_VEC_AOP_EXPR_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/vector/vec_expr.hpp>
#include <boost/numeric/mtl/operation/static_size.hpp>
#include <boost/numeric/mtl/operation/sfunctor.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/omp_size_type.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/unroll_size1.hpp>
#include <boost/numeric/mtl/utility/with_unroll1.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {

    namespace impl {

	template <unsigned long Index, unsigned long Max, typename SFunctor>
	struct assign
	{
	    typedef assign<Index+1, Max, SFunctor>     next;

	    template <typename E1, typename E2, typename Size>
	    static inline void apply(E1& first, const E2& second, Size i)
	    {
		SFunctor::apply( first(i+Index), second(i+Index) );
		next::apply( first, second, i );
	    }
	};

	template <unsigned long Max, typename SFunctor>
	struct assign<Max, Max, SFunctor>
	{
	    template <typename E1, typename E2, typename Size>
	    static inline void apply(E1& first, const E2& second, Size i)
	    {
		SFunctor::apply( first(i+Max), second(i+Max) );
	    }
	};
    }

// Generic assign operation expression template for vectors
// Model of VectorExpression
template <typename E1, typename E2, typename SFunctor>
struct vec_vec_aop_expr 
  :  public vec_expr< vec_vec_aop_expr<E1, E2, SFunctor> >
{
    typedef vec_expr< vec_vec_aop_expr<E1, E2, SFunctor> >  expr_base;
    typedef vec_vec_aop_expr<E1, E2, SFunctor>   self;
    typedef typename Collection<E1>::value_type  value_type;
    
    typedef typename Collection<E1>::size_type   size_type;
    typedef value_type reference_type ;

    typedef E1 first_argument_type ;
    typedef E2 second_argument_type ;
    
    vec_vec_aop_expr( first_argument_type& v1, second_argument_type const& v2, bool delay= false )
      : first(v1), second(v2), delayed_assign(delay)
    {
        bool compatible= mtl::vec::size(first) == mtl::vec::size(second) || (mtl::vec::size(first) == 0 && !traits::is_static<E1>::value);
        MTL_DEBUG_THROW_IF(!compatible,  incompatible_size());
	second.delay_assign();
    }
    
  private:
    void dynamic_assign(boost::mpl::false_) // Without unrolling
    {
	typedef typename traits::omp_size_type<size_type>::type size_type;
	size_type s= size_type(mtl::vec::size(first));

      #ifdef MTL_WITH_OPENMP
	# pragma omp parallel
	{
	    vampir_trace<8003> tracer;
	    #pragma omp for
	    for (size_type i= 0; i < s; ++i)
		SFunctor::apply( first(i), second(i) );
	}
      #else
	for (size_type i= 0; i < s; ++i)
	    SFunctor::apply( first(i), second(i) );
      #endif
    }

    void dynamic_assign(boost::mpl::true_) // With unrolling
    {
	typedef typename traits::omp_size_type<size_type>::type size_type;
	const size_type BSize= traits::unroll_size1<E1>::value0;
	size_type s= mtl::vec::size(first), sb= s / BSize * BSize;

      #ifdef MTL_WITH_OPENMP
	# pragma omp parallel
	{
	    vampir_trace<8003> tracer;
	    #pragma omp for
	    for (size_type i= 0; i < sb; i+= BSize)
		impl::assign<0, BSize-1, SFunctor>::apply(first, second, i);
	}
      #else
	for (size_type i= 0; i < sb; i+= BSize)
	    impl::assign<0, BSize-1, SFunctor>::apply(first, second, i);
      #endif

	for (size_type i= sb; i < s; i++) 
	    SFunctor::apply( first(i), second(i) );
    }    


    void assign(boost::mpl::false_)
    {
	vampir_trace<2017> tracer;
	// If target is constructed by default it takes size of source
	//int a= size(second);
	//int b= second;
	if (mtl::vec::size(first) == 0) 
            first.change_dim(mtl::size(second));

	// need to do more benchmarking before making unrolling default
	dynamic_assign(traits::with_unroll1<E1>());
    }

    void assign(boost::mpl::true_)
    {
	vampir_trace<1001> tracer;	
	// impl::assign<0, static_size<E1>::value-1, SFunctor>::apply(first, second); // Slower, at least on gcc
	for (size_type i= 0; i < mtl::vec::size(first); ++i) // Do an ordinary loop instead
	    SFunctor::apply(first(i), second(i));
    }

  public:
    ~vec_vec_aop_expr()
    {
	if (!delayed_assign)
	    assign(traits::is_static<E1>());
    }
    
    void delay_assign() const { delayed_assign= true; }

    template <typename EE1, typename EE2, typename SSFunctor>
    friend std::size_t size(const vec_vec_aop_expr<EE1, EE2, SSFunctor>& v);

    value_type& operator() (size_type i) const { 
	assert( delayed_assign );
	return SFunctor::apply(first(i), second(i));
    }

    value_type& operator[] (size_type i) const { return (*this)(i); }

    template <unsigned Offset>
    value_type& at(size_type i) const { 
	assert(delayed_assign);
	return SFunctor::apply(first(i+Offset), second(i+Offset));
    }
       
    first_argument_type const& first_argument() const { return first; }
    second_argument_type const& second_argument() const { return second; }

  private:
     first_argument_type&                first ;
     second_argument_type const&         second ;
     mutable bool                        delayed_assign;
} ; // vec_vec_aop_expr

template <typename E1, typename E2, typename SFunctor>
inline std::size_t size(const vec_vec_aop_expr<E1, E2, SFunctor>& v)
{
    MTL_DEBUG_THROW_IF( size(v.first) != 0 && size(v.first) != size(v.second), incompatible_size());
    return size(v.second);
}

} } // Namespace mtl::vector




#endif


