// Copyright 2006. Peter Gottschling, Matthias Troyer, Rolf Bonderer
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

#ifndef MATH_IDENTITY_INCLUDE
#define MATH_IDENTITY_INCLUDE

#include <boost/math/tools/config.hpp>
#include <limits>
#include <string>
#include <functional>

#include <boost/numeric/linear_algebra/operators.hpp>

namespace math {

template <typename Operation, typename Element>
struct identity_t {};

// TBD: Do we the case that the return type is different? Using std::unary_function?

// Additive identity of Element type is by default a converted 0
// However, for vectors one needs to know the dimension
// (and in parallel eventually also the distribution).
// Therefore, an element is passed as reference.
// It is strongly recommended to specialize this functor
// for better efficiency.
template <typename Element>
struct identity_t< add<Element>, Element > 
  : public std::binary_function<add<Element>, Element, Element>
{ 
    Element operator() (const add<Element>&, const Element& ref) const
    {
	Element tmp(ref);
	tmp= 0;
	return tmp;
    }
};

template <>
struct identity_t< add<std::string>, std::string > 
  : public std::binary_function<add<std::string>, std::string, std::string>
{ 
    std::string operator() (const add<std::string>&, const std::string&) const
    {
	return std::string();
    }
};

// Multiplicative identity of Element type is by default a converted 1
// Same comments as above.
// In contrast to additive identity, this default more likely to be wrong (e.g. matrices with all 1s)
template <typename Element>
struct identity_t< mult<Element>, Element > 
  : public std::binary_function<mult<Element>, Element, Element>
{ 
    Element operator() (const mult<Element>&, const Element& ref) const
    {
	Element tmp(ref);
	tmp= 1;
	return tmp;
    }
};


// Identity of max is minimal representable value, for standard types defined in numeric_limits
template <typename Element>
struct identity_t< max<Element>, Element > 
  : public std::binary_function<max<Element>, Element, Element>
{ 
    Element operator() (const max<Element>&, const Element& ) const
    {
	using std::numeric_limits;
	return numeric_limits<Element>::min();
    }
};

template <>
struct identity_t< max<float>, float > 
  : public std::binary_function<max<float>, float, float>
{ 
    float operator() (const max<float>&, const float& ) const
    {
	using std::numeric_limits;
	return -numeric_limits<float>::max();
    }
};

template <>
struct identity_t< max<double>, double > 
  : public std::binary_function<max<double>, double, double>
{ 
    double operator() (const max<double>&, const double& ) const
    {
	using std::numeric_limits;
	return -numeric_limits<double>::max();
    }
};


#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   
template <>
struct identity_t< max<long double>, long double > 
  : public std::binary_function<max<long double>, long double, long double>
{ 
    long double operator() (const max<long double>&, const long double& ) const
    {
	using std::numeric_limits;
	return -numeric_limits<long double>::max();
    }
};

#endif



// Identity of min is maximal representable value, for standard types defined in numeric_limits
template <typename Element>
struct identity_t< min<Element>, Element > 
  : public std::binary_function<min<Element>, Element, Element>
{ 
    Element operator() (const min<Element>&, const Element& ) const
    {
	using std::numeric_limits;
	return numeric_limits<Element>::max();
    }
};

// Identity of bit-wise and
template <typename Element>
struct identity_t< bitwise_and<Element>, Element > 
  : public std::binary_function<bitwise_and<Element>, Element, Element>
{ 
    Element operator() (const bitwise_and<Element>&, const Element&) const
    {
	return 0;
    }
};

// Identity of bit-wise or
template <typename Element>
struct identity_t< bitwise_or<Element>, Element > 
  : public std::binary_function<bitwise_or<Element>, Element, Element>
{ 
    Element operator() (const bitwise_or<Element>&, const Element&) const
    {
	return 0 - 1;
    }
};

#if 0 // ambiguous specialization
template <template <typename> class Operation, typename First, typename Second>
struct identity_t< Operation<std::pair<First, Second> >, std::pair<First, Second> >
{
    typedef std::pair<First, Second> pt;

    pt operator()(const Operation<pt>&, const pt& ref) const
    {
	return std::make_pair(identity(Operation<First>(), ref.first), identity(Operation<Second>(), ref.second));
    }
};
#endif

// Function is shorter than typetrait-like functor
template <typename Operation, typename Element>
inline Element identity(const Operation& op, const Element& v)
{
    return identity_t<Operation, Element>() (op, v);
}

#if 1
// I shouldn't do this (but as functor I'd need too many specializations)
template <template <typename> class Operation, typename First, typename Second>
inline std::pair<First, Second> identity(const Operation<std::pair<First, Second> >&, const std::pair<First, Second>& v)
{
    return std::pair<First, Second>(math::identity(Operation<First>(), v.first), math::identity(Operation<Second>(), v.second));
}
#endif

// Short-cut for additive identity
template <typename Element>
inline Element zero(const Element& v)
{
    return identity_t<math::add<Element>, Element>() (math::add<Element>(), v);
}


// Short-cut for multiplicative identity
template <typename Element>
inline Element one(const Element& v)
{
    return identity_t<math::mult<Element>, Element>() (math::mult<Element>(), v);
}


} // namespace math

#endif // MATH_IDENTITY_INCLUDE
