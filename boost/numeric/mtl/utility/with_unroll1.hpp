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

#ifndef MTL_TRAITS_WITH_UNROLL1_INCLUDE
#define MTL_TRAITS_WITH_UNROLL1_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Type trait for enabling one-dimensional unrolling, default is false.
template <typename Collection>
struct with_unroll1
  : boost::mpl::false_ {};

template <unsigned BSize, typename Vector>
struct with_unroll1<vec::unrolled1<BSize, Vector> >
  : boost::mpl::true_ {};


}} // namespace mtl::traits

#endif // MTL_TRAITS_WITH_UNROLL1_INCLUDE
