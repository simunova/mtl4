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

#ifndef MTL_FUSE_INCLUDE
#define MTL_FUSE_INCLUDE

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>

#include <boost/numeric/mtl/utility/is_lazy.hpp>
#include <boost/numeric/mtl/operation/fused_expr.hpp>

namespace mtl {

/// Fuse two lazy expressions (operator notation of fuse)
template <typename T, typename U>
typename boost::enable_if<boost::mpl::and_<traits::is_lazy<T>, traits::is_lazy<U> >, fused_expr<T, U> >::type
operator||(const T& x, const U& y)
{
    return fused_expr<T, U>(const_cast<T&>(x), const_cast<U&>(y));
}

/// Fuse two lazy expressions
template <typename T, typename U>
typename boost::enable_if<boost::mpl::and_<traits::is_lazy<T>, traits::is_lazy<U> >, fused_expr<T, U> >::type
fuse(const T& x, const U& y)
{
    return fused_expr<T, U>(const_cast<T&>(x), const_cast<U&>(y));
}



} // namespace mtl

#endif // MTL_FUSE_INCLUDE
