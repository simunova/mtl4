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

#ifndef MTL_IO_FUNCTOR_SYMBOL_INCLUDE
#define MTL_IO_FUNCTOR_SYMBOL_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace io {

template <typename Value> 
std::string functor_symbol(const sfunctor::conj<Value>&)
{  return "conj"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::imag<Value>&)
{  return "imag"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::real<Value>&)
{  return "real"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::identity<Value>&)
{  return "identity"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::abs<Value>&)
{  return "abs"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::sqrt<Value>&)
{  return "sqrt"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::square<Value>&)
{  return "square"; }

template <typename Value> 
std::string functor_symbol(const sfunctor::negate<Value>&)
{  return "-"; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::minus<Value1, Value2>&)
{  return "-"; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::plus<Value1, Value2>&)
{  return "+"; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::times<Value1, Value2>&)
{  return "*"; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::divide<Value1, Value2>&)
{  return "/"; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::assign<Value1, Value2>&)
{  return "="; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::plus_assign<Value1, Value2>&)
{  return "+="; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::minus_assign<Value1, Value2>&)
{  return "-="; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::times_assign<Value1, Value2>&)
{  return "*="; }

template <typename Value1, typename Value2> 
std::string functor_symbol(const sfunctor::divide_assign<Value1, Value2>&)
{  return "/="; }

}} // namespace mtl::io

#endif // MTL_IO_FUNCTOR_SYMBOL_INCLUDE
