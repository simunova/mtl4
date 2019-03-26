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

#ifndef MTL_PREDEFINED_MASKS_INCLUDE
#define MTL_PREDEFINED_MASKS_INCLUDE

#include <boost/numeric/mtl/recursion/bit_masking.hpp>

namespace mtl { namespace recursion {

    // Bitmasks: 
    const unsigned long morton_mask= generate_mask<true, 0, row_major, 0>::value,
			morton_z_mask= generate_mask<false, 0, row_major, 0>::value,
			doppled_2_row_mask= generate_mask<true, 1, row_major, 0>::value,
			doppled_2_col_mask= generate_mask<true, 1, col_major, 0>::value,
			doppled_4_row_mask= generate_mask<true, 2, row_major, 0>::value,
			doppled_4_col_mask= generate_mask<true, 2, col_major, 0>::value,
			doppled_16_row_mask= generate_mask<true, 4, row_major, 0>::value,
			doppled_16_col_mask= generate_mask<true, 4, col_major, 0>::value,
			doppled_z_16_row_mask= generate_mask<false, 4, row_major, 0>::value,
			doppled_z_16_col_mask= generate_mask<false, 4, col_major, 0>::value,
			doppled_32_row_mask= generate_mask<true, 5, row_major, 0>::value,
			doppled_32_col_mask= generate_mask<true, 5, col_major, 0>::value,
			doppled_z_32_row_mask= generate_mask<false, 5, row_major, 0>::value,
			doppled_z_32_col_mask= generate_mask<false, 5, col_major, 0>::value,
			doppled_64_row_mask= generate_mask<true, 6, row_major, 0>::value,
			doppled_64_col_mask= generate_mask<true, 6, col_major, 0>::value,
			doppled_z_64_row_mask= generate_mask<false, 6, row_major, 0>::value,
			doppled_z_64_col_mask= generate_mask<false, 6, col_major, 0>::value,
			doppled_128_row_mask= generate_mask<true, 7, row_major, 0>::value,
			doppled_128_col_mask= generate_mask<true, 7, col_major, 0>::value,
			doppled_z_128_row_mask= generate_mask<false, 7, row_major, 0>::value,
			doppled_z_128_col_mask= generate_mask<false, 7, col_major, 0>::value,
			shark_32_row_mask= generate_mask<true, 5, row_major, 1>::value,
			shark_32_col_mask= generate_mask<true, 5, col_major, 1>::value,
			shark_z_32_row_mask= generate_mask<false, 5, row_major, 1>::value,
			shark_z_32_col_mask= generate_mask<false, 5, col_major, 1>::value,
			shark_64_row_mask= generate_mask<true, 6, row_major, 1>::value,
			shark_64_col_mask= generate_mask<true, 6, col_major, 1>::value,
			shark_z_64_row_mask= generate_mask<false, 6, row_major, 1>::value,
			shark_z_64_col_mask= generate_mask<false, 6, col_major, 1>::value;


} // namespace recursion

// Export masks into MTL namespace

using recursion::morton_mask;
using recursion::morton_z_mask;
using recursion::doppled_2_row_mask;
using recursion::doppled_2_col_mask;
using recursion::doppled_4_row_mask;
using recursion::doppled_4_col_mask;
using recursion::doppled_16_row_mask;
using recursion::doppled_16_col_mask;
using recursion::doppled_z_16_row_mask;
using recursion::doppled_z_16_col_mask;
using recursion::doppled_32_row_mask;
using recursion::doppled_32_col_mask;
using recursion::doppled_z_32_row_mask;
using recursion::doppled_z_32_col_mask;
using recursion::doppled_64_row_mask;
using recursion::doppled_64_col_mask;
using recursion::doppled_z_64_row_mask;
using recursion::doppled_z_64_col_mask;
using recursion::doppled_128_row_mask;
using recursion::doppled_128_col_mask;
using recursion::doppled_z_128_row_mask;
using recursion::doppled_z_128_col_mask;
using recursion::shark_32_row_mask;
using recursion::shark_32_col_mask;
using recursion::shark_z_32_row_mask;
using recursion::shark_z_32_col_mask;
using recursion::shark_64_row_mask;
using recursion::shark_64_col_mask;
using recursion::shark_z_64_row_mask;
using recursion::shark_z_64_col_mask;

} // namespace mtl

#endif // MTL_PREDEFINED_MASKS_INCLUDE
