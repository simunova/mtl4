// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_TFUNCTOR_TFUNCTOR_MIXED_INCLUDE
#define MTL_TFUNCTOR_TFUNCTOR_MIXED_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <cmath>

namespace mtl { namespace tfunctor {

// Currently we define only scalar concepts
// Eventually this will be expanded like "scale"

/// Plus functor that stores left summand
template <typename Value1, typename Value2>
struct left_plus
{
    typedef typename Addable<Value1, Value2>::result_type result_type;	    
    explicit left_plus(const Value1& v1) : v1(v1) {}
    
    result_type operator() (const Value2& v2) const
    {
	return v1 + v2;
    }
private:
    Value1 v1; 
};

/// Plus functor that stores right summand
template <typename Value1, typename Value2>
struct right_plus
{
    typedef typename Addable<Value1, Value2>::result_type result_type;	    
    explicit right_plus(const Value2& v2) : v2(v2) {}
    
    result_type operator() (const Value1& v1) const
    {
	return v1 + v2;
    }
private:
    Value2 v2; 
};

/// Minus functor that stores left summand
template <typename Value1, typename Value2>
struct left_minus
{
    typedef typename Subtractable<Value1, Value2>::result_type result_type;	    
    explicit left_minus(const Value1& v1) : v1(v1) {}
    
    result_type operator() (const Value2& v2) const
    {
	return v1 - v2;
    }
private:
    Value1 v1; 
};

/// Minus functor that stores right summand
template <typename Value1, typename Value2>
struct right_minus
{
    typedef typename Subtractable<Value1, Value2>::result_type result_type;	    
    explicit right_minus(const Value2& v2) : v2(v2) {}
    
    result_type operator() (const Value1& v1) const
    {
	return v1 - v2;
    }
private:
    Value2 v2; 
};

/// Minimum functor that stores left operand
/** Result type is right operand. 
    This way the min of a scalar and a vector keeps the vector type. **/
template <typename Value1, typename Value2>
struct left_min
{
    typedef Value2 result_type;	    
    explicit left_min(const Value1& v1) : v1(v1) {}
    
    result_type operator() (const Value2& v2) const
    {
	return v1 < v2 ? Value2(v1) : v2;
    }
private:
    Value1 v1; 
};

/// Minimum functor that stores right operand
/** Result type is left operand. 
    This way the min of a vector and a scalar keeps the vector type. **/
template <typename Value1, typename Value2>
struct right_min
{
    typedef Value1 result_type;	    
    explicit right_min(const Value2& v2) : v2(v2) {}
    
    result_type operator() (const Value1& v1) const
    {
	return v1 < v2 ? v1 : Value1(v2);
    }
private:
    Value2 v2; 
};

/// Maximum functor that stores left operand
/** Result type is right operand. 
    This way the max of a scalar and a vector keeps the vector type. **/
template <typename Value1, typename Value2>
struct left_max
{
    typedef Value2 result_type;	    
    explicit left_max(const Value1& v1) : v1(v1) {}
    
    result_type operator() (const Value2& v2) const
    {
	return v1 < v2 ? v2 : Value2(v1);
    }
private:
    Value1 v1; 
};

/// Maximum functor that stores right operand
/** Result type is left operand. 
    This way the max of a vector and a scalar keeps the vector type. **/
template <typename Value1, typename Value2>
struct right_max
{
    typedef Value1 result_type;	    
    explicit right_max(const Value2& v2) : v2(v2) {}
    
    result_type operator() (const Value1& v1) const
    {
	return v1 < v2 ? Value1(v2) : v1;
    }
private:
    Value2 v2; 
};

/// Raise Value1 to the power of Value2.
/** Value2 is bound in the functor and Value1 passes to the operator. **/
template <typename Value1, typename Value2>
struct pow_by
{
    typedef Value1 result_type; // not perfect but works for all C++98 overloads
    
    explicit pow_by(const Value2& exp) : exp(exp) {}
    
    result_type operator() (const Value1& base) const
    {
        using std::pow;
        return pow(base, exp);
    }

    Value2 value() const { return exp; }
private:
    Value2 exp; 
};


}} // namespace mtl::tfunctor

#endif // MTL_TFUNCTOR_TFUNCTOR_MIXED_INCLUDE
