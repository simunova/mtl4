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

#ifndef MTL_TRAITS_LINEAR_OPERATOR_INCLUDE
#define MTL_TRAITS_LINEAR_OPERATOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Type trait returning type of linear operator projecting from vector space Vector1 to vector space Vector2
/** The operator is dense. **/
template <typename Vector1, typename Vector2>
struct linear_operator {};	

template <typename Value1, typename Para1, typename Value2, typename Para2> 
struct linear_operator<mtl::dense_vector<Value1, Para1>, mtl::dense_vector<Value2, Para2> >
{
    typedef mtl::mat::dense2D<Value1>   type;
};


}} // namespace mtl::traits

#endif // MTL_TRAITS_LINEAR_OPERATOR_INCLUDE
