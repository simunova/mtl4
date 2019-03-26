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

#ifndef MTL_STATIC_SIZE_INCLUDE
#define MTL_STATIC_SIZE_INCLUDE

#include <boost/numeric/mtl/operation/static_num_rows.hpp>
#include <boost/numeric/mtl/operation/static_num_cols.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

/// Number of rows times columns given at compile time
/** General declaration, relies on static_num_rows and static_num_cols. **/
template <typename Collection>
struct static_size {
	vampir_trace<1007> tracer;
    typedef typename static_num_rows<Collection>::type type; // Might not always work
    static const type value= static_num_rows<Collection>::value * static_num_cols<Collection>::value;
};

} // namespace mtl

#endif // MTL_STATIC_SIZE_INCLUDE
