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

#ifndef MTL_SFUNCTOR_INCLUDE
#define MTL_SFUNCTOR_INCLUDE

#include <cmath>
#include <complex>

#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/concept/static_functor.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/type_traits.hpp>

namespace mtl { namespace sfunctor {

template <typename Value1, typename Value2>
struct plus
{
    typedef const Value1&                                 first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef typename Addable<Value1, Value2>::result_type result_type;

    static inline result_type apply(const Value1& v1, const Value2& v2)
    {
  return v1 + v2;
    }

    result_type operator() (const Value1& v1, const Value2& v2) const
    {
  vampir_trace<23> tracer;
  return v1 + v2;
    }
};
    
template <typename Value1, typename Value2>
struct minus
{
    typedef const Value1&                                 first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef typename Subtractable<Value1, Value2>::result_type result_type;

    static inline result_type apply(const Value1& v1, const Value2& v2)
    {
  return v1 - v2;
    }

    result_type operator() (const Value1& v1, const Value2& v2) const
    {
  vampir_trace<24> tracer;
  return v1 - v2;
    }
};

template <typename Value1, typename Value2>
struct times
{
    typedef const Value1&                                 first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef typename Multiplicable<Value1, Value2>::result_type result_type;

    static inline result_type apply(const Value1& v1, const Value2& v2)
    {
  return v1 * v2;
    }

    result_type operator() (const Value1& v1, const Value2& v2) const
    {
  vampir_trace<25> tracer;
  return v1 * v2;
    }
};

template <typename Value1, typename Value2>
struct divide
{
    typedef const Value1&                                 first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef typename Divisible<Value1, Value2>::result_type result_type;

    static inline result_type apply(const Value1& v1, const Value2& v2)
    {
  return v1 / v2;
    }

    result_type operator() (const Value1& v1, const Value2& v2) const
    {
  vampir_trace<26> tracer;
  return v1 / v2;
    }
};

template <typename Value1, typename Value2>
struct assign
{
    typedef Value1&                                       first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef Value1&                                       result_type;

    static inline result_type apply(Value1& v1, const Value2& v2)
    {
  return v1= Value1(v2);
    }

    result_type operator() (Value1& v1, const Value2& v2) const
    {
  vampir_trace<27> tracer;
  return v1= v2;
    }
};
    
template <typename Value1, typename Value2>
struct plus_assign
{
    typedef Value1&                                       first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef Value1&                                       result_type;

    static inline result_type apply(Value1& v1, const Value2& v2)
    {
  return v1+= v2;
    }

    result_type operator() (Value1& v1, const Value2& v2) const
    {
  vampir_trace<28> tracer;
  return v1+= v2;
    }
};
    
template <typename Value1, typename Value2>
struct minus_assign
{
    typedef Value1&                                       first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef Value1&                                       result_type;

    static inline result_type apply(Value1& v1, const Value2& v2)
    {
  return v1-= v2;
    }

    result_type operator() (Value1& v1, const Value2& v2) const
    {
  vampir_trace<29> tracer;
  return v1-= v2;
    }
};

template <typename Value1, typename Value2>
struct times_assign
{
    typedef Value1&                                       first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef Value1&                                       result_type;

    static inline result_type apply(Value1& v1, const Value2& v2)
    {
  return v1*= v2;
    }

    result_type operator() (Value1& v1, const Value2& v2) const
    {
  vampir_trace<30> tracer;
  return v1*= v2;
    }
};

template <typename Value1, typename Value2>
struct divide_assign
{
    typedef Value1&                                       first_argument_type;
    typedef const Value2&                                 second_argument_type;
    typedef Value1&                                       result_type;

    static inline result_type apply(Value1& v1, const Value2& v2)
    {
  return v1/= v2;
    }

    result_type operator() (Value1& v1, const Value2& v2) const
    {
  vampir_trace<31> tracer;
  return v1/= v2;
    }
};


// Might be helpful for surplus functor arguments
template <typename Value>
struct identity
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v)
    {
  return v;
    }

    result_type operator() (const Value& v) const
    {
  vampir_trace<32> tracer;
  return v;
    }
};


template <typename Value>
struct abs
{
    typedef const Value&                                  argument_type;
    typedef typename Magnitude<Value>::type               result_type;

    static inline result_type apply(const Value& v)
    {            
  using std::abs;
  return abs(v);
    }

    result_type operator() (const Value& v)  const
    {
  vampir_trace<33> tracer; 
  return apply(v); 
    }
};

template <typename Value>
struct sqrt
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v)
    {            
  using std::sqrt;
  return sqrt(v);
    }

    result_type operator() (const Value& v)  const
    {
  vampir_trace<34> tracer;
  return apply(v); 
    }
};

template <typename Value>
struct square
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v)
    {            
  return v * v;
    }

    result_type operator() (const Value& v) const
    {
  vampir_trace<35> tracer;
  return apply(v);
    }
};


template <typename Value>
struct negate
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) { return -v;  }
    result_type operator() (const Value& v) const 
    {
  vampir_trace<36> tracer;
  return -v;
    }
};

template <typename Value>
struct acos
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::acos;
        return acos(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

# ifdef MTL_WITH_MATH_ELEVEN
template <typename Value>
struct acosh
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::acosh;
        return acosh(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};
# endif

template <typename Value>
struct asin
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::asin;
        return asin(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

# ifdef MTL_WITH_MATH_ELEVEN
template <typename Value>
struct asinh
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::asinh;
        return asinh(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};
# endif

template <typename Value>
struct atan
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::atan;
        return atan(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

# ifdef MTL_WITH_MATH_ELEVEN
template <typename Value>
struct atanh
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::atanh;
        return atanh(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};
# endif


template <typename Value>
struct cos
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::cos;
        return cos(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

template <typename Value>
struct cosh
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::cosh;
        return cosh(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

template <typename Value>
struct sin
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::sin;
        return sin(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

template <typename Value>
struct sinh
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::sinh;
        return sinh(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

template <typename Value>
struct tan
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::tan;
        return tan(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

template <typename Value>
struct tanh
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::tanh;
        return tanh(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

// Rounding functions

template <typename Value>
struct ceil
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    {
        return apply(v, boost::is_integral<Value>());
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }

private:
    static inline result_type apply(const Value& v, boost::integral_constant<bool, false>)
    {
        using std::ceil;
        return ceil(v);
    }
    
    // return value directly for integer values
    static inline result_type apply(const Value& v, boost::integral_constant<bool, true>)
    {
        return v;
    };
};

template <typename Value>
struct floor
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    {
        return apply(v, boost::is_integral<Value>());
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }

private:
    static inline result_type apply(const Value& v, boost::integral_constant<bool, false>)
    {
        using std::floor;
        return floor(v);
    }
    
    // return value directly for integer values
    static inline result_type apply(const Value& v, boost::integral_constant<bool, true>)
    {
        return v;
    };
};

# ifdef MTL_WITH_MATH_ELEVEN    

    template <typename Value>
    struct round
    {
        typedef const Value&                                  argument_type;
        typedef Value                                         result_type;

        static inline result_type apply(const Value& v) 
        {
            return apply(v, boost::is_integral<Value>());
        }
        result_type operator() (const Value& v) const 
        {
            return apply(v);
        }

    private:
        static inline result_type apply(const Value& v, boost::integral_constant<bool, false>)
        {
            using std::round;
            return round(v);
        }
        
        // return value directly for integer values
        static inline result_type apply(const Value& v, boost::integral_constant<bool, true>)
        {
            return v;
        };
    };

    template <typename Value>
    struct trunc
    {
        typedef const Value&                                  argument_type;
        typedef Value                                         result_type;

        static inline result_type apply(const Value& v) 
        {
            return apply(v, boost::is_integral<Value>());
        }
        result_type operator() (const Value& v) const 
        {
            return apply(v);
        }

    private:
        static inline result_type apply(const Value& v, boost::integral_constant<bool, false>)
        {
            using std::trunc;
            return trunc(v);
        }
        
        // return value directly for integer values
        static inline result_type apply(const Value& v, boost::integral_constant<bool, true>)
        {
            return v;
        };
    };

# endif

// Logarithmic functions

template <typename Value>
struct log
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::log;
        return log(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Value>
struct log2
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::log2;
        return log2(v);  
    }
    result_type operator() (const Value& v) const 
    {
        return apply(v);
    }
};
# endif

template <typename Value>
struct log10
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::log10;
        return log10(v);  
    }
    result_type operator() (const Value& v) const 
    {
        return apply(v);
    }
};

// Exponential functions

template <typename Value>
struct exp
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::exp;
        return exp(v);  
    }
    result_type operator() (const Value& v) const 
    {
        return apply(v);
    }
};

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Value>
struct exp2
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::exp2;
        return exp2(v);  
    }
    result_type operator() (const Value& v) const 
    {
        return apply(v);
    }
};
#endif

template <typename Value>
struct exp10
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::exp;
        return exp(v * 2.302585092994045684017991454684364207601101488628772976033);  
    }
    result_type operator() (const Value& v) const 
    {
        return apply(v);
    }
};

// Inverse square root functions

template <typename Value>
struct rsqrt
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::pow;
        return pow(v, -0.5);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

// Error functions

# ifdef MTL_WITH_MATH_ELEVEN    
template <typename Value>
struct erf
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::erf;
        return erf(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};

template <typename Value>
struct erfc
{
    typedef const Value&                                  argument_type;
    typedef Value                                         result_type;

    static inline result_type apply(const Value& v) 
    { 
        using std::erfc;
        return erfc(v);  
    }
    result_type operator() (const Value& v) const 
    {
  return apply(v);
    }
};
# endif

/// Compose functors \p F and \p G, i.e. compute f(g(x)).
/** Functors must be models of StaticUnaryFunctor,
    StaticUnaryFunctor<G>::result_type must be convertible to
    StaticUnaryFunctor<F>::argument_type.
    Under these conditions compose<F, G> will be a model of StaticUnaryFunctor.
**/
template <typename F, typename G>
struct compose
{
    typedef typename StaticUnaryFunctor<G>::argument_type argument_type;
    typedef typename StaticUnaryFunctor<F>::result_type   result_type;
    
    static inline result_type apply(argument_type x)
    {
  return F::apply(G::apply(x));
    }

    result_type operator()(argument_type x) 
    {
  vampir_trace<37> tracer;
  return apply(x);
    }
};


/// Compose functors \p F and \p G with G in F's first argument, i.e. compute f(g(x), y).
/** F/G must be models of StaticBinaryFunctor/StaticUnaryFunctor,
    StaticUnaryFunctor<G>::result_type must be convertible to
    StaticBinaryFunctor<F>::first_argument_type.
    Under these conditions compose_first<F, G> will be a model of StaticBinaryFunctor.
**/
template <typename F, typename G>
struct compose_first
{
    typedef typename StaticUnaryFunctor<G>::argument_type         first_argument_type;
    typedef typename StaticBinaryFunctor<F>::second_argument_type second_argument_type;
    typedef typename StaticBinaryFunctor<F>::result_type          result_type;
    
    static inline result_type apply(first_argument_type x, second_argument_type y)
    {
  return F::apply(G::apply(x), y);
    }

    result_type operator()(first_argument_type x, second_argument_type y)
    {
  vampir_trace<38> tracer;
  return apply(x, y);
    }
};


/// Compose functors \p F and \p G with G in F's second argument, i.e. compute f(x, g(y)).
/** F/G must be models of StaticBinaryFunctor/StaticUnaryFunctor,
    StaticUnaryFunctor<G>::result_type must be convertible to
    StaticBinaryFunctor<F>::second_argument_type.
    Under these conditions compose_second<F, G> will be a model of StaticBinaryFunctor.
**/
template <typename F, typename G>
struct compose_second
{
    typedef typename StaticBinaryFunctor<F>::first_argument_type  first_argument_type;
    typedef typename StaticUnaryFunctor<G>::argument_type         second_argument_type;
    typedef typename StaticBinaryFunctor<F>::result_type          result_type;
    
    static inline result_type apply(first_argument_type x, second_argument_type y)
    {
  return F::apply(x, G::apply(y));
    }

    result_type operator()(first_argument_type x, second_argument_type y)
    {
  vampir_trace<39> tracer;
  return apply(x, y);
    }
};

/// Compose functors \p F, \p G, and \p H with G/H in F's first/second argument, i.e. compute f(g(x), h(y)).
/** F/G must be models of StaticBinaryFunctor/StaticUnaryFunctor,
    StaticUnaryFunctor<G>::result_type must be convertible to
    StaticBinaryFunctor<F>::first_argument_type and
    StaticUnaryFunctor<H>::result_type must be convertible to
    StaticBinaryFunctor<F>::second_argument_type.
    Under these conditions compose_both<F, G, H> will be a model of StaticBinaryFunctor.
**/
template <typename F, typename G, typename H>
struct compose_both
{
    typedef typename StaticUnaryFunctor<G>::argument_type         first_argument_type;
    typedef typename StaticUnaryFunctor<H>::argument_type         second_argument_type;
    typedef typename StaticBinaryFunctor<F>::result_type          result_type;
    
    static inline result_type apply(first_argument_type x, second_argument_type y)
    {
  return F::apply(G::apply(x), H::apply(y));
    }

    result_type operator()(first_argument_type x, second_argument_type y)
    {
  vampir_trace<40> tracer;
  return apply(x, y);
    }
};

/// Compose unary functor \p F with binary functor \p G, i.e. compute f(g(x, y)).
/** F/G must be models of StaticUnaryFunctor/StaticBinaryFunctor,
    StaticBinaryFunctor<G>::result_type must be convertible to
    StaticUnaryFunctor<F>::argument_type.
    Under these conditions compose_binary<F, G> will be a model of StaticBinaryFunctor.
**/
template <typename F, typename G>
struct compose_binary
{
    typedef typename StaticBinaryFunctor<G>::first_argument_type  first_argument_type;
    typedef typename StaticBinaryFunctor<G>::second_argument_type second_argument_type;
    typedef typename StaticUnaryFunctor<F>::result_type           result_type;
    
    static inline result_type apply(first_argument_type x, second_argument_type y)
    {
  return F::apply(G::apply(x, y));
    }

    result_type operator()(first_argument_type x, second_argument_type y)
    {
  vampir_trace<41> tracer;
  return apply(x, y);
    }
};


/// Templatized example of composition, computes l_2 norm in 2D, i.e. sqrt(abs(x*x + y*y))
template <typename T>
struct l_2_2D
  : public compose_binary<sqrt<typename abs<T>::result_type>, 
        compose_binary<abs<T>, 
           compose_both<plus<T, T>, 
                  square<T>, 
                  square<T>  > 
                                        > 
                         >
{};

}} // namespace mtl::sfunctor

#endif // MTL_SFUNCTOR_INCLUDE
