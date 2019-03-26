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

#ifndef MTL_DIMENSIONS_INCLUDE
#define MTL_DIMENSIONS_INCLUDE

#include <iostream>
#include <cassert>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

namespace mtl {

// dimension is a type for declaring matrix dimensions 
// num_rows() and num_cols() return the number or rows and columns
// is_static says whether it is declared at compile time or not

// Compile time version
namespace fixed
{

    /// Compile-time dimensions
    template <std::size_t Rows, std::size_t Cols>
    struct dimensions
    {
	typedef std::size_t  size_type;

        static size_type const Num_Rows= Rows;
        static size_type const Num_Cols= Cols;

	// To have the same interface as fixed
#ifndef NDEBUG
	/// Constructor does not need arguments but if given they are compared against the template arguments in debug mode
	explicit dimensions(size_type r= Rows, size_type c= Cols) 
	{
	    assert(r == Rows); assert(c == Cols); 
	}
#else
	explicit dimensions(size_type, size_type) {}
	explicit dimensions(size_type) {}
	explicit dimensions() {}
#endif

	size_type num_rows() const { return Rows; } ///< Number of rows
	size_type num_cols() const { return Cols; } ///< Number of columns

	/// To check whether dimensions are static
	static bool const is_static= true;

	/// Transposed dimension (type)
	typedef dimensions<Cols, Rows> transposed_type;
	transposed_type transpose() const 
	{ 
	    return transposed_type(); 
	}
    };

    /// Output of dimensions
    template <std::size_t R, std::size_t C>
    inline std::ostream& operator<< (std::ostream& stream, dimensions<R, C>) 
    {
	return stream << R << 'x' << C; 
    }

} // namespace fixed

namespace non_fixed
{
    /// Run-time dimensions
    struct dimensions
    {
	typedef std::size_t  size_type;

	/// Constructor 
	dimensions(size_type r= 0, size_type c= 0) : r(r), c(c) {}
	
	/// Assignment
	dimensions& operator=(const dimensions& x) 
	{
	    r= x.r; c= x.c; return *this; 
	}
	size_type num_rows() const { return r; } ///< Number of rows
	size_type num_cols() const { return c; } ///< Number of columns

	/// Transposed dimension
	typedef dimensions transposed_type;
	transposed_type transpose() 
	{ 
	    return transposed_type(c, r); 
	}

	/// To check whether dimensions are static
	static bool const is_static= false;
    protected:
	size_type r, c;
    };

    /// Output of dimensions
    inline std::ostream& operator<< (std::ostream& stream, dimensions d) 
    {
	return stream << d.num_rows() << 'x' << d.num_cols(); 
    }

} // namespace non_fixed

} // namespace mtl

#endif // MTL_DIMENSIONS_INCLUDE
