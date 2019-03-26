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

#ifndef MTL_IRANGE_INCLUDE
#define MTL_IRANGE_INCLUDE

#include <limits>
#include <iostream>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/detail/base_cursor.hpp>


namespace mtl { 

    /// Maximal index
    const std::size_t imax= std::numeric_limits<std::size_t>::max();

    /// Class to define a half open index ranges 
    class irange
    {
      public:

        typedef std::size_t                    size_type;
	typedef detail::base_cursor<size_type> cursor_type;

        /// Create an index range of [start, finish)
        explicit irange(size_type start, size_type finish) : my_start(start), my_finish(finish) {}

        /// Create an index range of [0, finish)
        explicit irange(size_type finish) : my_start(0), my_finish(finish) {}

        /// Create an index range of [0, imax), i.e. all indices
        irange() : my_start(0), my_finish(imax) {}

        /// Set the index range to [start, finish)
	irange& set(size_type start, size_type finish) 
	{
	    my_start= start; my_finish= finish; return *this;
	}

        /// Set the index range of [0, finish)
	irange& set(size_type finish) 
	{
	    my_start= 0; my_finish= finish; return *this;
	}

        /// Decrease finish, i.e. [start, finish) -> [start, finish-1)
	irange& operator--() 
	{
	    --my_finish; return *this;
	}
	
        /// First index in range
        size_type start() const { return my_start; } 
        /// Past-end index in range
        size_type finish() const { return my_finish; }
        /// Number of indices
        size_type size() const { return my_finish > my_start ? my_finish - my_start : 0; }

	/// Whether the range is empty
        bool empty() const { return my_finish <= my_start; }

        /// First index in range
	cursor_type begin() const { return my_start; } 
        /// Past-end index in range
        cursor_type end() const { return my_finish; }

	/// Maps integers [0, size()) to [start(), finish())
	/** Checks index in debug mode. Inverse of from_range. **/
	size_type to_range(size_type i) const
	{
	    MTL_DEBUG_THROW_IF(is_negative(i) || i >= size(), index_out_of_range());
	    return my_start + i;
	}

	/// Maps integers [start(), finish()) to [0, size())
	/** Checks index in debug mode. **/
	size_type from_range(size_type i) const
	{
	    MTL_DEBUG_THROW_IF(i < my_start || i >= my_finish, index_out_of_range());
	    return i - my_start;
	}

	/// Wether index is in range
	bool in_range(size_type i) const
	{
	    return i >= my_start && i < my_finish;
	}

	friend std::ostream& operator<<(std::ostream& os, const irange& ir)
	{   return os << "mtl::irange(" << ir.my_start << ", " << ir.my_finish << ")"; }

      private:
        size_type my_start, my_finish;
    };
    /// All index in range
    namespace {
	// problems with multiple cpp projects, which include mtl multiple times
	static const irange iall(0, imax);
	//inline irange iall() { return irange(); }
    }

    /// Intersection of two ranges
    irange inline intersection(irange const& r1, irange const& r2)
    {
	return irange(std::max(r1.start(), r2.start()), std::min(r1.finish(), r2.finish()));
    }


} // namespace mtl



#endif // MTL_IRANGE_INCLUDE
