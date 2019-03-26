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

#ifndef MTL_ENABLE_IF_INCLUDE
#define MTL_ENABLE_IF_INCLUDE

#include <boost/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>

namespace mtl { namespace traits {

template <typename Value, typename Type = void>
struct enable_if_matrix
  : boost::enable_if<is_matrix<Value>, Type>
{};

template <typename Value, typename Type = void>
struct enable_if_vector
  : boost::enable_if<is_vector<Value>, Type>
{};

template <typename Value, typename Type = void>
struct enable_if_scalar
  : boost::enable_if<is_scalar<Value>, Type>
{};

template <typename Value, typename Type>
struct lazy_enable_if_matrix
  : boost::lazy_enable_if<is_matrix<Value>, Type>
{};

template <typename Value, typename Type>
struct lazy_enable_if_vector
  : boost::lazy_enable_if<is_vector<Value>, Type>
{};

template <typename Value, typename Type>
struct lazy_enable_if_scalar
  : boost::lazy_enable_if<is_scalar<Value>, Type>
{};

}} // namespace mtl

#endif // MTL_ENABLE_IF_INCLUDE
