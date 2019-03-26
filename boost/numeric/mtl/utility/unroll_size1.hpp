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

#ifndef MTL_TRAITS_UNROLL_SIZE1_INCLUDE
#define MTL_TRAITS_UNROLL_SIZE1_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Type trait for one-dimensional unrolling, default is 4.
template <typename Collection>
struct unroll_size1
{
    static const unsigned value0= 4;
};

template <unsigned BSize, typename Vector>
struct unroll_size1<vec::unrolled1<BSize, Vector> >
{
    static const unsigned value0= BSize;
};


}} // namespace mtl::traits

#endif // MTL_TRAITS_UNROLL_SIZE1_INCLUDE
