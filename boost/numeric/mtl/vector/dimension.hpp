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

#ifndef MTL_DIMENSION_INCLUDE
#define MTL_DIMENSION_INCLUDE

namespace mtl { namespace vec {

// Compile time version
namespace fixed {

    /// Compile-time dimension
    template <std::size_t Size>
    struct dimension
    {
	typedef std::size_t  size_type;
	static size_type const value= Size;
	friend inline size_type size(const dimension&) { return dimension::value; } ///< Size

	/// To check whether it is static
	static bool const is_static= true;
    };
}

namespace non_fixed {

    /// Run-time dimension
    struct dimension
    {
	typedef std::size_t  size_type;
	
	static size_type const value= 0; // for compatibility
	dimension(size_type v= 0) : my_size(v) {} /// Constructor with default zero
	friend inline size_type size(const dimension& d) { return d.my_size; } ///< Size

	/// To check whether it is static
	static bool const is_static= false;
      protected:
	size_type my_size;
    };
}

}} // namespace mtl::vec

#endif // MTL_DIMENSION_INCLUDE
