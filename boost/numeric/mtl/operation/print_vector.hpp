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

#ifndef MTL_PRINT_VECTOR_INCLUDE
#define MTL_PRINT_VECTOR_INCLUDE

#include <iostream>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>

namespace mtl { namespace vec {

template <typename Vector>
std::ostream& print_vector(Vector const& vector, std::ostream& out= std::cout, int width= 0, int precision= 0)
{
    using mtl::vec::size;
    out << '{' << mtl::vec::size(vector) 
	<< (traits::is_row_major< typename OrientedCollection<Vector>::orientation >::value ? "R" : "C") 
	<< "}[" ;
    for (size_t r = 0; r < mtl::vec::size(vector); ++r) {
	out.fill (' '); 
	if (width) out.width (width); 
	// out.flags (std::ios_base::right);
	if (precision) out.precision(precision); 
	out << vector[r] << (r+1 < mtl::vec::size(vector) ? "," : "");
    }
    return out << ']';
}

}} // namespace mtl::vector

#endif // MTL_PRINT_VECTOR_INCLUDE
