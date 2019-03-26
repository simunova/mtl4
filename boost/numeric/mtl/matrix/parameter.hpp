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

#ifndef MTL_MATRIX_PARAMETERS_INCLUDE
#define MTL_MATRIX_PARAMETERS_INCLUDE

#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/detail/index.hpp>
#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>

namespace mtl { namespace mat {

/// Type for bundling template parameters of common matrix types
/** OnStack = true can only be used with fixed::dimensions.
    \sa \ref matrix_parameters
    \sa \ref tuning_fsize
    \sa \ref tuning_sizetype **/
template <typename Orientation= row_major, 
	  typename Index= index::c_index,
	  typename Dimensions= mtl::non_fixed::dimensions,
	  bool OnStack= mtl::traits::is_static<Dimensions>::value,
	  typename SizeType= std::size_t>
struct parameters 
{
    typedef Orientation orientation;
    typedef Index       index;
    typedef Dimensions  dimensions;
    static bool const   on_stack= OnStack;
    typedef SizeType    size_type;

    // Matrix dimensions must be known at compile time to be on the stack
    // MTL_STATIC_ASSERT(( !on_stack || dimensions::is_static ), "Types to be stored on stack must provide static size.");
};

/// Short-cut to define parameters with unsigned and defaults otherwise
typedef parameters<row_major, index::c_index, mtl::non_fixed::dimensions, false, unsigned> unsigned_parameters;

}} // namespace mtl::matrix

#endif // MTL_MATRIX_PARAMETERS_INCLUDE
