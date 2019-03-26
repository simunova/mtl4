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

#ifndef MTL_FUSED_EXPR_INCLUDE
#define MTL_FUSED_EXPR_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/numeric/mtl/operation/index_evaluator.hpp>
#include <boost/numeric/mtl/utility/index_evaluatable.hpp>
#include <boost/numeric/mtl/utility/assert.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

/// Expression template for fusing other expression 
template <typename T, typename U>
struct fused_expr
{
    typedef boost::mpl::and_<traits::forward_index_evaluatable<T>, traits::forward_index_evaluatable<U> >   forward;
    typedef boost::mpl::and_<traits::backward_index_evaluatable<T>, traits::backward_index_evaluatable<U> > backward;

    fused_expr(T& first, U& second) 
      : first(first), second(second), delayed_assign(false) 
    {
	first.delay_assign(); second.delay_assign();
    }
 
    ~fused_expr() { if (!delayed_assign) eval(forward(), backward()); }

    void delay_assign() const { delayed_assign= true; }    

    template <typename TT, typename UU>
    void forward_eval_loop(const TT& const_first_eval, const UU& const_second_eval, boost::mpl::false_)
    {	
	vampir_trace<6001> tracer;
	// hope there is a more elegant way; copying the arguments causes errors due to double destructor evaluation
	TT& first_eval= const_cast<TT&>(const_first_eval);  
	UU& second_eval= const_cast<UU&>(const_second_eval);
	MTL_CRASH_IF(mtl::size(first_eval) != mtl::size(second_eval), "Incompatible size!");	

	for (std::size_t i= 0, s= size(first_eval); i < s; i++) {
	    first_eval(i); second_eval(i);
	}	
    }

    template <typename TT, typename UU>
    void forward_eval_loop(const TT& const_first_eval, const UU& const_second_eval, boost::mpl::true_)
    {	
	vampir_trace<6002> tracer;
	// hope there is a more elegant way; copying the arguments causes errors due to double destructor evaluation
	TT& first_eval= const_cast<TT&>(const_first_eval);  
	UU& second_eval= const_cast<UU&>(const_second_eval);
	MTL_CRASH_IF(mtl::vec::size(first_eval) != mtl::vec::size(second_eval), "Incompatible size!");	

	const std::size_t s= size(first_eval), sb= s >> 2 << 2;

	for (std::size_t i= 0; i < sb; i+= 4) {
	    first_eval.template at<0>(i); second_eval.template at<0>(i);
	    first_eval.template at<1>(i); second_eval.template at<1>(i);
	    first_eval.template at<2>(i); second_eval.template at<2>(i);
	    first_eval.template at<3>(i); second_eval.template at<3>(i);
	}

	for (std::size_t i= sb; i < s; i++) {
	    first_eval(i); second_eval(i);
	}
    }

    // Forward evaluation dominates backward
    template <bool B2>
    void eval(boost::mpl::true_, boost::mpl::bool_<B2>)
    {
#ifdef MTL_LAZY_LOOP_WO_UNROLL
	typedef boost::mpl::false_                                                                              to_unroll;
#else
	typedef boost::mpl::and_<traits::unrolled_index_evaluatable<T>, traits::unrolled_index_evaluatable<U> > to_unroll;
#endif
	// Currently lazy evaluation is only available on vector expressions, might change in the future
	// std::cout << "Forward evaluation\n";
	forward_eval_loop(index_evaluator(first), index_evaluator(second), to_unroll()); 
    }

    template <typename TT, typename UU>
    void backward_eval_loop(const TT& const_first_eval, const UU& const_second_eval, boost::mpl::false_)
    {	
	vampir_trace<6003> tracer;
	// hope there is a more elegant way; copying the arguments causes errors due to double destructor evaluation
	TT& first_eval= const_cast<TT&>(const_first_eval);  
	UU& second_eval= const_cast<UU&>(const_second_eval);
	MTL_CRASH_IF(mtl::size(first_eval) != mtl::size(second_eval), "Incompatible size!");	

	for (std::size_t i= size(first_eval); i-- > 0; ) {
	    // std::cout << "i is " << i << "\n";
	    first_eval(i); second_eval(i);
	}
    }

    template <typename TT, typename UU>
    void backward_eval_loop(const TT& const_first_eval, const UU& const_second_eval, boost::mpl::true_)
    {	
	vampir_trace<6004> tracer;
	// hope there is a more elegant way; copying the arguments causes errors due to double destructor evaluation
	TT& first_eval= const_cast<TT&>(const_first_eval);  
	UU& second_eval= const_cast<UU&>(const_second_eval);
	MTL_CRASH_IF(mtl::size(first_eval) != mtl::size(second_eval), "Incompatible size!");	

	std::size_t s= size(first_eval), i= s-1, m= s % 4;
	for (; m; m--) {
	    // std::cout << "i is " << i << "\n";
	    first_eval(i); second_eval(i--);
	}
	for(long j= i - 3; j >= 0; j-= 4) {
	    // std::cout << "i is " << j+3 << ".." << j << "\n";
	    first_eval.template at<3>(j); second_eval.template at<3>(j);
	    first_eval.template at<2>(j); second_eval.template at<2>(j);
	    first_eval.template at<1>(j); second_eval.template at<1>(j);
	    first_eval.template at<0>(j); second_eval.template at<0>(j);
	}
    }

    // Backward evaluation if forward isn't possible
    void eval(boost::mpl::false_, boost::mpl::true_)
    {
	// std::cout << "Backward evaluation\n";
	backward_eval_loop(index_evaluator(first), index_evaluator(second), boost::mpl::true_());
    }

    // Sequential evaluation
    void eval(boost::mpl::false_, boost::mpl::false_)
    { 
	// std::cout << "Non-fused evaluation\n";
	evaluate_lazy(first); evaluate_lazy(second); 
    }

    T& first;
    U& second;
    mutable bool                        delayed_assign;
};


} // namespace mtl

#endif // MTL_FUSED_EXPR_INCLUDE
