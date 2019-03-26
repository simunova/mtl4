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

#ifndef MTL_STRING_TO_ENUM_INCLUDE
#define MTL_STRING_TO_ENUM_INCLUDE

#include <string>
#include <boost/numeric/mtl/operation/size.hpp>

namespace mtl {

/** Searches string \p s in list \p l of strings and returns enum
    
    List \p l is given as array of const char*, which is the easiest to
    initialize.  The search is case sensitive, thus (de)-capitalize
    your string upfront, e.g., with boost::to_lower().
    If the string is not found then a runtime_error is thrown.
**/
template <typename EnumType, typename Array>
EnumType inline string_to_enum(const std::string& s, const Array& l, EnumType)
{
    std::size_t i;
    for (i= 0; i < size(l) && std::string(l[i]) != s; i++) {}
    MTL_THROW_IF(i == size(l), runtime_error("Search string not found"));
    return EnumType(i);
}

} // namespace mtl

#endif // MTL_STRING_TO_ENUM_INCLUDE
