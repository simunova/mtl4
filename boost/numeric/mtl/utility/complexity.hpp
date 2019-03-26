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

#ifndef MTL_COMPLEXITY_INCLUDE
#define MTL_COMPLEXITY_INCLUDE

#include <boost/mpl/int.hpp>
#include <boost/mpl/max.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/less_equal.hpp>

// This file contains types to characterize time complexities of operations or traversals
// Different traversals over given collections have different run time
// Algorithms implementable with different traversals can use these complexities to dispatch

// Complexities in a order, which of course will not be changed
// The underlying MPL definitions might be modified to add finer grained distinctions



namespace mtl { 
/// Namespace to define complexity classes
namespace complexity_classes {

/// Constant complexity where the value is expected to be in cache with high probability.
/** Special type for traversals to distinguish between strided or random memory access with 'constant' 
    (but slow) memory access and consecutive memory access with a good change that only one element
    per cache line must be load from memory
**/
struct cached : boost::mpl::int_<1> {};

/// Constant complexity 
struct constant : boost::mpl::int_<2> {};

/// Logarithmic complexity 
struct log_n : boost::mpl::int_<4> {};

/// Polynomial logarithm, i.e. log^k n
struct polylog_n : boost::mpl::int_<5> {};

/// Product of linear and cached
struct linear_cached : boost::mpl::int_<21> {};

/// Linear
struct linear : boost::mpl::int_<22> {};

/// Logarithm times linear, i.e. n * log n
struct n_log_n : boost::mpl::int_<24> {};

/// Polynomial logarithm times linear, i.e. n * log^k n
struct n_polylog_n : boost::mpl::int_<25> {};

struct quadratic : boost::mpl::int_<41> {};

/// All complexities larger than quadratic (< infinite) including n^2 log^k n
struct polynomial : boost::mpl::int_<200> {};

/// Infinite time complexity, which usually means that the operation or traversal is not available
struct infinite : boost::mpl::int_<1000> {};

/// Adding complexities of two operations is the maximal complexity of both operations
template <typename X, typename Y> 
struct plus : boost::mpl::if_< boost::mpl::less<X, Y>, Y, X> {};

/// Minimal complexity class
template <typename X, typename Y> 
struct min : boost::mpl::if_< boost::mpl::less<X, Y>, X, Y> {};


namespace detail
{
    // specializations on first argument

    // polynomial is the most frequent result, otherwise explicit definition later 
    template <typename X, typename Y> struct times 
    {
	typedef polynomial type;
    };

    template <typename Y> struct times<cached, Y> 
    {
	typedef Y type; 
    };

    template <typename Y> struct times<constant, Y> 
    {
	typedef Y type; 
    };   

    template <> struct times<log_n, log_n> 
    {
	typedef polylog_n type; 
    };   
    
    template <> struct times<log_n, polylog_n> : times<log_n, log_n> {};

    template <> struct times<log_n, linear_cached> 
    {
	typedef n_log_n type; 
    };   
    
    template <> struct times<log_n, linear> : times<log_n, linear_cached> {};
    
    template <> struct times<log_n, n_log_n> 
    {
	typedef n_polylog_n type; 
    };   
    
    template <> struct times<log_n, n_polylog_n> : times<log_n, n_log_n> {};
    
    template <> struct times<polylog_n, polylog_n> 
    {
	typedef polylog_n type; 
    };   
    
    template <> struct times<polylog_n, linear_cached> 
    {
	typedef n_polylog_n type; 
    };   

    template <> struct times<polylog_n, linear> : times<polylog_n, linear_cached> {};
    
    template <> struct times<polylog_n, n_log_n> : times<polylog_n, linear_cached> {};
    
    template <> struct times<polylog_n, n_polylog_n> : times<polylog_n, linear_cached> {};

    template <> struct times<linear_cached, linear_cached> 
    {
	typedef quadratic type; 
    };   
    
    template <> struct times<linear_cached, linear> : times<linear_cached, linear_cached> {};

    template <> struct times<linear, linear> : times<linear_cached, linear_cached> {};

} // namespace detail

/// Product of complexities
// Multiplication needs to be defined explicitly
// At least is symmetric, so we only consider X <= Y
template <typename X, typename Y> 
struct times
    : boost::mpl::if_< 
            boost::mpl::less<X, Y>
          , detail::times<X, Y>
          , detail::times<Y, X>
    >::type
{};

// Specializations on second argument (if were ordered)
// Done here to avoid ambiguities

template <typename X> 
struct times<X, infinite>
{
    typedef infinite type;
};

template <typename X> struct times<infinite, X> : times<X, infinite> {};

}} // namespace mtl::complexity

#define MTL_PRINT_COMPLEXITY(TYPE, STRING) \
inline std::ostream& operator<< (std::ostream& os, mtl::complexity_classes::TYPE) \
{                                                 \
    return os << STRING;                          \
}

MTL_PRINT_COMPLEXITY(cached, "cached constant complexity")
MTL_PRINT_COMPLEXITY(constant, "constant complexity")
MTL_PRINT_COMPLEXITY(log_n, "logarithmic complexity")
MTL_PRINT_COMPLEXITY(polylog_n, "poly-logarithmic complexity")
MTL_PRINT_COMPLEXITY(linear_cached, "cached linear complexity")
MTL_PRINT_COMPLEXITY(linear, "linear complexity")
MTL_PRINT_COMPLEXITY(n_log_n, "n log n complexity")
MTL_PRINT_COMPLEXITY(n_polylog_n, "n poly-log n complexity")
MTL_PRINT_COMPLEXITY(quadratic, "quadratic complexity")
MTL_PRINT_COMPLEXITY(polynomial, "polynomial complexity")
MTL_PRINT_COMPLEXITY(infinite, "infinite complexity")



#endif // MTL_COMPLEXITY_INCLUDE
