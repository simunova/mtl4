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

#ifndef MTL_REDUCTION_FUNCTORS_INCLUDE
#define MTL_REDUCTION_FUNCTORS_INCLUDE

#include <cmath>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/operation/squared_abs.hpp>

namespace mtl { namespace vec {

struct one_norm_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::zero;
	value= zero(value);
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	using std::abs;
	value+= abs(x);
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	value+= value2;
    }

    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};


// sub-optimal if abs is not needed
struct two_norm_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::zero;
	value= zero(value);
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	using mtl::squared_abs;
	value+= squared_abs(x);
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	value+= value2;
    }

    // After reduction compute square root
    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	using std::sqrt;
	return sqrt(value);
    }
};

// same as two-norm without the root at the end
struct unary_dot_functor
  : two_norm_functor
{
    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};

struct infinity_norm_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::zero;
	value= zero(value);
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	using std::abs; using std::max;
	value= max(value, Value(abs(x)));
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	using std::abs; using std::max;
	value= max(value, Value(abs(value2)));
    }

    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};


struct sum_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::zero;
	value= zero(value);
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	value+= x;
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	value+= value2;
    }

    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};


struct product_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::one;
	value= one(value);
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	value*= x;
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	value*= value2;
    }

    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};


struct max_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::identity; 
	value= math::identity(math::max<Value>(), value); // ADL doesn't work here in g++ 4.4
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	value= math::max<Value>()(value, x);
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	value= math::max<Value>()(value, value2);
    }

    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};


struct min_functor
{
    template <typename Value>
    static inline void init(Value& value)
    {
	using math::identity; 
	value= math::identity(math::min<Value>(), value); // ADL doesn't work here in g++ 4.4
    }

    template <typename Value, typename Element>
    static inline void update(Value& value, const Element& x)
    {    
	value= math::min<Value>()(value, x);
    }

    template <typename Value>
    static inline void finish(Value& value, const Value& value2)
    {
	value= math::min<Value>()(value, value2);
    }

    template <typename Value>
    static inline Value post_reduction(const Value& value)
    {
	return value;
    }
};


}} // namespace mtl::vector

#endif // MTL_REDUCTION_FUNCTORS_INCLUDE
