// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_SRANGE_INCLUDE
#define MTL_SRANGE_INCLUDE

#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>

namespace mtl {

class srange
{
  public:
    typedef std::size_t size_type;

    /// Create a strided index range of [start, finish)
    explicit srange(size_type start, size_type finish, size_type stride) 
      : my_start(start), my_finish(finish), my_stride(stride) {}
    
    /// First index in range
    size_type start() const { return my_start; } 
    /// Past-end index in range
    size_type finish() const { return my_finish; }
    /// Stride
    size_type stride() const { return my_stride; }

    /// Number of indices
    size_type size() const { return my_finish > my_start ? my_finish - my_start / my_stride : 0; }

    /// Maps integers [0, size()) to [start(), start()+stride(), ..., finish())
    /** Checks index in debug mode. Inverse of from_range. **/
    size_type to_range(size_type i) const
    {
	MTL_DEBUG_THROW_IF(is_negative(i) || i >= size(), index_out_of_range());
	return my_start + i * my_stride;
    }

    /// Maps integers [start(), finish()) to [0, size())
    /** Checks index in debug mode. **/
    size_type from_range(size_type i) const
    {
	MTL_DEBUG_THROW_IF(i < my_start || i >= my_finish, index_out_of_range());
	MTL_DEBUG_THROW_IF((i - my_start) % my_stride != 0, runtime_error("Index not on stride."));
	return i - my_start;
    }

    /// Wether index is in range
    bool in_range(size_type i) const
    {
	return i >= my_start && i < my_finish && (i - my_start) % my_stride == 0;
    }

   private:
    size_type my_start, my_finish, my_stride;
};

} // namespace mtl

#endif // MTL_SRANGE_INCLUDE
