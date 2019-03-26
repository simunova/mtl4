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

#ifndef MTL_SHIFT_BLOCKS_INCLUDE
#define MTL_SHIFT_BLOCKS_INCLUDE

#include <boost/numeric/mtl/operation/shift_block_detail.hpp>

namespace mtl { namespace operations {

// Shift blocks in an 1D array to remove unnecessary holes
// inserting holes in other places where needed (e.g. for inserting new values)
//
// Block 'i' is the half-open interval [starts[i], ends[i]) in data
// It will be copied into [new_starts[i], ...) in place
// Blocks are ordered: start[i] <= start[i+1]
// Data between blocks are considered holes and can be overwritten
//
template <typename Size, typename Starts, typename NewStarts, typename Ends, typename Data>
void shift_blocks(Size blocks, Starts const& starts, NewStarts const& new_starts, 
		  Ends const& ends, Data& data)
{
    for (Size i = 0; i < blocks; ) {
	detail::copy_blocks_forward(i, blocks, starts, new_starts, ends, data);
	detail::copy_blocks_backward(i, blocks, starts, new_starts, ends, data);
    }
}


}} // namespace mtl::operations

#endif // MTL_SHIFT_BLOCKS_INCLUDE
