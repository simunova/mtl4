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


#ifndef MATH_ACCUMULATE_INCLUDE
#define MATH_ACCUMULATE_INCLUDE

#include <concepts>

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/new_concepts.hpp>

namespace math {

// Dispatching between simple and unrolled version
template <std::InputIterator Iter, std::CopyConstructible Value, typename Op>
     requires std::Callable2<Op, Value, Value>
           && std::CopyAssignable<Value, std::Callable2<Op, Value, Value>::result_type>
Value inline accumulate(Iter first, Iter last, Value init, Op op)
{
    // std::cout << "Simple accumulate\n";
    for (; first != last; ++first)
	init= op(init, Value(*first));
    return init;
}


template <std::RandomAccessIterator Iter, std::CopyConstructible Value, typename Op>
    requires std::Callable2<Op, Value, Value>
          && std::CopyAssignable<Value, std::Callable2<Op, Value, Value>::result_type>
          && Commutative<Op, Value> 
          && Monoid<Op, Value> 
          && std::Convertible<Monoid<Op, Value>::identity_result_type, Value>
Value inline accumulate(Iter first, Iter last, Value init, Op op)
{
    // std::cout << "Unrolled accumulate\n";
    typedef typename std::RandomAccessIterator<Iter>::difference_type difference_type;
    Value            t0= identity(op, init), t1= identity(op, init), 
	             t2= identity(op, init), t3= init;
    difference_type  size= last - first, bsize= size >> 2 << 2, i;
    
    for (i= 0; i < bsize; i+= 4) {
	t0= op(t0, Value(first[i]));
	t1= op(t1, Value(first[i+1]));
	t2= op(t2, Value(first[i+2]));
	t3= op(t3, Value(first[i+3]));
    }
    for (; i < size; i++)
	t0= op(t0, Value(first[i]));

    t0= op(t0, t1), t2= op(t2, t3), t0= op(t0, t2);
    return t0;
}

} // namespace math

#endif // MATH_ACCUMULATE_INCLUDE
