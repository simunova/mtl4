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

#ifndef MTL_ISET_INCLUDE
#define MTL_ISET_INCLUDE

#include <vector>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/numeric/mtl/utility/push_back_comma_inserter.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl {

    /// Class for arbitrary index sets 
    class iset
    {
      public:
	/// Size type
        typedef std::size_t size_type;

	void check(size_type MTL_DEBUG_ARG(i)) const
	{ MTL_DEBUG_THROW_IF(is_negative(i) || i >= indices.size(), index_out_of_range()); }

	iset() {} ///< Default constructor
	
	/// Construct from array
	template <long Size>
	iset(const size_type (&array)[Size]) : indices(&array[0], &array[Size]) {}

	/// Constructor from std::vector; size_type must be identic
	explicit iset(const std::vector<size_type>& src) : indices(src) {}

#     ifdef MTL_WITH_MOVE
	/// Constructor from std::vector; size_type must be identic
	explicit iset(std::vector<size_type>&& src) : indices(std::move(src)) {}
#    endif

	/// Assign comma-separated list
	template <typename Source>
	typename boost::enable_if<boost::is_integral<Source>, push_back_comma_inserter<iset> >::type
	operator=(const Source& src)
	{
	    indices.clear();
	    indices.push_back(src);
	    return push_back_comma_inserter<iset>(*this);
	}

	/// Return i-th index
	size_type operator[](size_type i) const { check(i); return indices[i]; }

	/// Insert a new index at the end
	void push_back(size_type i) { indices.push_back(i); }

	/// Size of the set
	size_type size() const { return indices.size(); }

	/// Equality
	bool operator==(const iset& that) const
	{
	    if (size() != that.size()) 
		return false;
	    for (size_type i= 0; i < size(); i++)
		if (indices[i] != that.indices[i]) 
		    return false;
	    return true;
	}

	bool operator!=(const iset& that) const { return !(*this == that); } ///< Inequality

	/// Print iset
	friend std::ostream& operator<<(std::ostream& os, const iset& is)
	{   os << "{";
	    for (size_type i= 0; i < is.indices.size(); i++) {
		os << is.indices[i];
		if (i+1 < is.indices.size()) 
		    os << ", ";
	    }
	    return os << "}"; 
	}

      private:
	std::vector<size_type> indices;
    };

    /// Size of an index set (as free function)
    std::size_t inline size(const iset& i) 
    { return i.size(); }

} // namespace mtl

#endif // MTL_ISET_INCLUDE
