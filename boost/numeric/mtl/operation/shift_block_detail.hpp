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

#ifndef MTL_SHIFT_BLOCKS_DETAIL_INCLUDE
#define MTL_SHIFT_BLOCKS_DETAIL_INCLUDE

#include <algorithm>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace operations { namespace detail {

// &v[i] is replaced by &v[0]+i to enable past-end addresses for STL copy
// otherwise MSVC complaines

template <typename Size, typename Starts, typename NewStarts, typename Ends, typename Data>
inline void copy_blocks_forward(Size& i, Size blocks, Starts const& starts, NewStarts const& new_starts, 
			 Ends const& ends, Data& data)
{
    using std::copy;
	vampir_trace<1005> tracer;
    // Copy forward as long as blocks are not shifted 
    for (; i < blocks && starts[i] >= new_starts[i]; ++i) 	
	if (starts[i] > new_starts[i])
	    copy(&data[0] + starts[i], &data[0] + ends[i], &data[0] + new_starts[i]);
}

template <typename Size, typename Starts, typename NewStarts, typename Ends, typename Data>
inline void copy_blocks_backward(Size& i, Size blocks, Starts const& starts, NewStarts const& new_starts, 
			  Ends const& ends, Data& data)
{
    using std::copy;
    using std::copy_backward;
	vampir_trace<1006> tracer;
    Size first = i;
    // find first block to be copied forward (or end)
    while (i < blocks && starts[i] < new_starts[i]) ++i;

    for (Size j = i; j-- > first; )
	if (ends[j] <= new_starts[j])
	    copy(&data[0] + starts[j], &data[0] + ends[j], &data[0] + new_starts[j]);
	else
	    copy_backward(&data[0] + starts[j], &data[0] + ends[j], &data[0] + (new_starts[j]+ends[j]-starts[j]));
}

}}} // namespace mtl::operations::detail

#endif // MTL_SHIFT_BLOCKS_DETAIL_INCLUDE
