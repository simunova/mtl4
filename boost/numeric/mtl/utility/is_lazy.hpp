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

#ifndef MTL_TRAITS_IS_LAZY_INCLUDE
#define MTL_TRAITS_IS_LAZY_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/mpl/and.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Type trait for determining wether \p T is to be evaluated lazily
template <typename T>
struct is_lazy : boost::mpl::false_ {};

template <typename T, typename U, typename Assign>
struct is_lazy<lazy_assign<T, U, Assign> > : boost::mpl::true_ {};

template <typename T, typename U>
struct is_lazy<fused_expr<T, U> > 
  : boost::mpl::and_<is_lazy<T>, is_lazy<U> > 
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_IS_LAZY_INCLUDE
