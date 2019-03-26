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

#ifndef MTL_VECTOR_PARAMETERS_INCLUDE
#define MTL_VECTOR_PARAMETERS_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/vector/dimension.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>

namespace mtl { namespace vec {

/// This type exist only for bundling template parameters (to reduce typing)
/** OnStack = true can only be used with fixed::dimension.
    \sa \ref tuning_fsize
    \sa \ref tuning_sizetype **/
template <typename Orientation= col_major, 
	  typename Dimension= non_fixed::dimension,
	  bool OnStack= mtl::traits::is_static<Dimension>::value,
	  typename SizeType= std::size_t>
struct parameters 
{
    typedef Orientation orientation;
    typedef Dimension   dimension;
    static const bool   on_stack= OnStack;
    typedef SizeType    size_type;

    // Vector dimension must be known at compile time to be on the stack
    MTL_STATIC_ASSERT(( !on_stack || dimension::is_static ), "Types to be stored on stack must provide static size.");
};


}} // namespace mtl::vector

#endif // MTL_VECTOR_PARAMETERS_INCLUDE
