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

#ifndef MTL_STD_OUTPUT_OPERATOR_INCLUDE
#define MTL_STD_OUTPUT_OPERATOR_INCLUDE

#ifdef MTL_HAS_STD_OUTPUT_OPERATOR

#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace std {

    /// Print standard vector
    /** Only available when compiled with macro flag MTL_HAS_STD_OUTPUT_OPERATOR
	to avoid (reduce) conflicts with other libraries. **/
    template <typename T, typename Alloc>
    inline ostream& operator<< (ostream& os, vector<T, Alloc> const&  v)
    {
	os << '[';
	for (size_t r = 0; r < v.size(); ++r)
	    os << v[r] << (r < v.size() - 1 ? "," : "");
	return os << ']';
    }

    /// Print standard map
    /** Only available when compiled with macro flag MTL_HAS_STD_OUTPUT_OPERATOR
	to avoid (reduce) conflicts with other libraries. **/
    template <typename Key, typename Data, typename Compare, typename Alloc>
    inline ostream& operator<< (ostream& os, map<Key, Data, Compare, Alloc> const&  m)
    {
	if (m.empty()) return os << "{}";
	typedef typename map<Key, Data, Compare, Alloc>::const_iterator iter_type;
	os << '{';
	for (iter_type it= m.begin(), end= m.end(); it != end; ++it)
	    os << it->first << ": " << it->second << ", ";
	return os << "\b\b} ";
    }

    /// Print standard pair
    /** Only available when compiled with macro flag MTL_HAS_STD_OUTPUT_OPERATOR
	to avoid (reduce) conflicts with other libraries. **/
    template <typename T, typename U>
    inline ostream& operator<< (ostream& os, pair<T, U> const& p)
    {
	return os << '(' << p.first << ',' << p.second << ')';
    }

} // namespace std

#endif // MTL_HAS_STD_OUTPUT_OPERATOR

#endif // MTL_STD_OUTPUT_OPERATOR_INCLUDE
