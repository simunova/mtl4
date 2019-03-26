#ifndef algebraic_functions_include
#define algebraic_functions_include

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <boost/config/concept_macros.hpp> 

#ifdef __GXX_CONCEPTS__
#  include <concepts>
#  include <bits/iterator_concepts.h>
#endif

#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>

// Pure algebraic functions (mainly for applying concepts)

namespace mtl {


// {Op, Element} must be a Magma
// The result of the operation must be EqualityComparable
// Closed2EqualityComparable requires the comparability of any combination of
// Element and result_type.
// It is thus more limiting than necessary, well slightly.
template <typename Op, typename Element>
#if 0
  _GLIBCXX_WHERE( math::Magma<Op, Element>
	    && std::EqualityComparable< math::Magma<Op, Element>::result_type > )
#endif
  _GLIBCXX_WHERE( math::Closed2EqualityComparable<Op, Element> && math::Magma<Op, Element> )
inline bool equal_results(const Element& v1a, const Element& v1b, 
			  const Element& v2a, const Element& v2b, Op op) 
{
#if 0
    // If we use EqualityComparable for Element we cannot compare the results
    // directly since Magma only assumes that the results are convertible to Element
    // which does not imply that they are comparable.
    // Assigning to temporaries will the conversion and the temporaries are comparable.

    Element res1= op(v1a, v1b), res2= op(v2a, v2b);
    return res1 == res2;
#endif

    return op(v1a, v1b) == op(v2a, v2b);
}


// Same for AdditiveMagma
template <typename Element>
_GLIBCXX_WHERE ( math::AdditiveMagma<Element> 
	   // && math::Closed2EqualityComparable< math::add<Element>, Element > )
	   && std::EqualityComparable< math::AdditiveMagma<Element>::addition_result_type > )
inline bool equal_add_results(const Element& v1a, const Element& v1b, 
			      const Element& v2a, const Element& v2b) 
{
    return v1a + v1b == v2a + v2b;
}




// {Op, Element} must be a Monoid
// The result of the operation must be EqualityComparable
template <typename Op, typename Element>
  _GLIBCXX_WHERE( math::Closed2EqualityComparable<Op, Element> && math::Monoid<Op, Element> )
inline bool identity_pair(const Element& v1, const Element& v2, Op op) 
{
    using math::identity;
    return op(v1, v2) == identity(op, v1) ;
}


// {Op, Element} must be a SemiGroup
template <typename Op, typename Element, typename Exponent>
  _GLIBCXX_WHERE( math::SemiGroup<Op, Element> 
            && std::Integral<Exponent> )  
Element recursive_multiply_and_square(const Element& a, Exponent n, Op op) 
{
    if (n <= 0) throw "In recursive_multiply_and_square: exponent must greater than 0";

    Exponent half= n >> 1;

    // If half is 0 then n must be 1 and the result is a
    if (half == 0)
	return a;

    // compute power of downward rounded exponent and square the result
    Element value= recursive_multiply_and_square(a, half, op);
    value= op(value, value);

    // if odd another multiplication with a is needed
    if (n & 1) 
	value= op(value, a);
    return value;
} 


// {Op, Element} must be a SemiGroup
template <typename Op, typename Element, typename Exponent>
  _GLIBCXX_WHERE( math::SemiGroup<Op, Element> 
            && std::Integral<Exponent> )  
inline Element multiply_and_square_horner(const Element& a, Exponent n, Op op) 
{
    if (n <= 0) throw "In multiply_and_square_horner: exponent must be greater than 0";

    // Set mask to highest bit
    Exponent mask= 1 << (8 * sizeof(mask) - 1);

    // If this is a negative number right shift can insert 1s instead of 0s -> infinite loop
    // Therefore we take the 2nd-highest bit
    if (mask < 0)
	mask= 1 << (8 * sizeof(mask) - 2);

    // find highest 1 bit
    while(!(bool)(mask & n)) mask>>= 1;

    Element value= a;
    for (mask>>= 1; mask; mask>>= 1) {
	value= op(value, value);
	if (n & mask) 
	    value= op(value, a);
    }
    return value;
}


// {Op, Element} must be a Monoid
template <typename Op, typename Element, typename Exponent>
  _GLIBCXX_WHERE( math::Monoid<Op, Element> 
            && std::Integral<Exponent> ) 
inline Element multiply_and_square(const Element& a, Exponent n, Op op) 
{
    // Same as the simpler form except that the first multiplication is made before 
    // the loop and one squaring is saved this way
    if (n < 0) throw "In multiply_and_square: negative exponent";

    using math::identity;
    Element value= bool(n & 1) ? Element(a) : Element(identity(op, a)), square= a;

    for (n>>= 1; n > 0; n>>= 1) {
	square= op(square, square); 
	if (n & 1) 
	    value= op(value, square);
    }
    return value;  
} 


// {Op, Element} must be a Monoid
template <typename Op, typename Element, typename Exponent>
  _GLIBCXX_WHERE( math::Monoid<Op, Element> 
            && std::Integral<Exponent> ) 
inline Element multiply_and_square_simple(const Element& a, Exponent n, Op op) 
{
    if (n < 0) throw "In multiply_and_square: negative exponent";

    using math::identity;
    Element value= identity(op, a), square= a;
    for (; n > 0; n>>= 1) {
	if (n & 1) 
	    value= op(value, square);
	square= op(square, square); 
    }
    return value;  
} 


template <typename Iter, typename Value, typename Op>
  _GLIBCXX_WHERE( std::ForwardIterator<Iter> 
                  && std::Convertible<Value, std::ForwardIterator<Iter>::value_type>
		  && math::Magma<Op, std::ForwardIterator<Iter>::value_type> )
typename std::ForwardIterator<Iter>::value_type 
inline accumulate_simple(Iter first, Iter last, Value init, Op op)
{
    typedef typename std::RandomAccessIterator<Iter>::value_type value_type;
    value_type        t0= init;
    
    // std::cout << "accumulate_simple\n";
    for (; first != last; ++first)
	t0= op(t0, *first);
    return t0;
}

template <typename Iter, typename Value, typename Op>
  _GLIBCXX_WHERE( std::RandomAccessIterator<Iter> 
	    && std::Convertible<Value, std::RandomAccessIterator<Iter>::value_type>
	    && math::CommutativeMonoid<Op, std::RandomAccessIterator<Iter>::value_type> )
typename std::RandomAccessIterator<Iter>::value_type 
inline accumulate_unrolled(Iter first, Iter last, Value init, Op op)
{
    typedef typename std::RandomAccessIterator<Iter>::value_type value_type;
    typedef typename std::RandomAccessIterator<Iter>::difference_type difference_type;
    value_type        t0= identity(op, init), t1= identity(op, init), t2= identity(op, init), t3= init;
    difference_type size= last - first, bsize= size >> 2 << 2, i;
    
    for (i= 0; i < bsize; i+= 4) {
	t0= op(t0, first[i]);
	t1= op(t1, first[i+1]);
	t2= op(t2, first[i+2]);
	t3= op(t3, first[i+3]);
    }
    for (; i < size; i++)
	t0= op(t0, first[i]);
    return op(op(t0, t1), op(t2, t3));
}



// Dispatching between simple and unrolled version
template <typename Iter, typename Value, typename Op>
  _GLIBCXX_WHERE( std::ForwardIterator<Iter> 
                  && std::Convertible<Value, std::ForwardIterator<Iter>::value_type>
		  && math::Magma<Op, std::ForwardIterator<Iter>::value_type> )
typename std::ForwardIterator<Iter>::value_type 
inline accumulate(Iter first, Iter last, Value init, Op op)
{
    std::cout << "Simple accumulate\n";
    return accumulate_simple(first, last, init, op);
}


template <typename Iter, typename Value, typename Op>
  _GLIBCXX_WHERE( std::RandomAccessIterator<Iter> 
	          && std::Convertible<Value, std::RandomAccessIterator<Iter>::value_type>
		  && math::CommutativeMonoid<Op, std::RandomAccessIterator<Iter>::value_type> )
typename std::RandomAccessIterator<Iter>::value_type 
inline accumulate(Iter first, Iter last, Value init, Op op)
{
    std::cout << "Unrolled accumulate\n";
    return accumulate_unrolled(first, last, init, op);
}




// {Op, Element} must be a PartiallyInvertibleMonoid
// Element and results must be EqualityComparable
// Only 
template <typename Op, typename Element>
  _GLIBCXX_WHERE( math::Closed2EqualityComparable<Op, Element> 
	    && math::PartiallyInvertibleMonoid<Op, Element> )
inline int algebraic_division(const Element& v1, const Element& v2, Op op)
{
    using math::identity; using math::inverse; using math::is_invertible;

    if (!is_invertible(op, v2))
	throw "In algebraic division: v2 must be invertible!\n";

    // Temporaries to avoid redundant operations
    Element id= identity(op, v1),     // Identity
            iv2= inverse(op, v2),     // Inverse of v2
   	    tmp(v1);                             // Copy of v1, will be lessened until < id

    int counter= 0;
    for (; tmp != id; counter++) 
	tmp= op(tmp, iv2);
    return counter;
}


// {Op, Element} must be a PartiallyInvertibleMonoid
// Element must be LessThanComparable
// Under construction w.r.t. semantic requirements, introduction of ordered group needed
template <typename Op, typename Element>
  _GLIBCXX_WHERE( math::Closed2LessThanComparable<Op, Element> 
	    && math::PartiallyInvertibleMonoid<Op, Element> )
inline int ordered_algebraic_division(const Element& v1, const Element& v2, Op op)
{
    using math::identity; using math::inverse; using math::is_invertible;

    if (!is_invertible(op, v2))
	throw "In algebraic division: v2 must be invertible!\n";

    // Temporaries to avoid redundant operations
    Element id= identity(op, v1),     // Identity
            iv2= inverse(op, v2),     // Inverse of v2
   	    tmp(v1);                             // Copy of v1, will be lessened until < id
    
    if (v1 <= id) return 0;
    int counter= 0;
    for (; tmp > id; counter++) 
	tmp= op(tmp, iv2);
    // counter only correct if tmp == id
    if (tmp < id) 
	counter--;
    return counter;
}




#if 0

// {Iter*, Op} must be a CommutativeMonoid
struct sortedAccumulate_t {
  template <class Iter, class Op, class Comp>
  typename enable_if<glas::is_associative<typename std::iterator_traits<Iter>::value_type, Op>::value 
                       && glas::is_commutative<typename std::iterator_traits<Iter>::value_type, Op>::value, 
		     typename std::iterator_traits<Iter>::value_type>::type
  operator() (Iter first, Iter last, Op op, Comp comp) {
    std::cout << "sortedAccumulate_t\n";
    typedef typename std::iterator_traits<Iter>::value_type value_type;
    std::vector<value_type> tmp(first, last);
    std::sort(tmp.begin(), tmp.end(), comp); 
    return std::accumulate(tmp.begin(), tmp.end(), glas::identity<value_type, Op>()(), op); }
} sortedAccumulate;

// {Iter*, Op} must be a Monoid
struct unsortedAccumulate_t {
  template <class Iter, class Op>
  typename std::iterator_traits<Iter>::value_type
  operator() (Iter first, Iter last, Op op) {
    std::cout << "unsortedAccumulate_t\n";
    typedef typename std::iterator_traits<Iter>::value_type value_type;
    return std::accumulate(first, last, glas::identity<value_type, Op>()(), op); }

  // Only for Compability
  template <class Iter, class Op, class Comp>
  typename std::iterator_traits<Iter>::value_type
  operator() (Iter first, Iter last, Op op, Comp) {
    return operator() (first, last, op); }
} unsortedAccumulate;
    
// {Iter*, Op} must be a Monoid
template <class Iter, class Op, class Comp>
inline typename std::iterator_traits<Iter>::value_type
trySortedAccumulate(Iter first, Iter last, Op op, Comp comp) {
  typedef typename std::iterator_traits<Iter>::value_type value_type;
  typename if_type<glas::is_associative<value_type, Op>::value && glas::is_commutative<value_type, Op>::value,  
                   sortedAccumulate_t, unsortedAccumulate_t>::type  accumulate;
  // alternatively checking for structure
  //   typename if_type<glas::is_commutative_monoid<value_type, Op>::value,  
  //                    sortedAccumulate_t, unsortedAccumulate_t>::type  accumulate;
  return accumulate(first, last, op, comp);
}

#endif // 0

} // namespace mtl

#endif // algebraic_functions_include
