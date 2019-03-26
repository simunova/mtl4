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

#ifndef MTL_TRAITS_DOMAIN_INCLUDE
#define MTL_TRAITS_DOMAIN_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>

namespace mtl { namespace traits {

// Helper
template <typename Matrix, bool IsMatrix>
struct domain_aux {};

template <typename Matrix>
struct domain_aux<Matrix, true>
{
    typedef mtl::vec::dense_vector<typename Collection<Matrix>::value_type, vec::parameters<> > type;
};

/// Type trait returning domain type of linear operator (matrix)
template <typename Matrix>
struct domain
  : domain_aux<Matrix, is_matrix<Matrix>::value>
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_DOMAIN_INCLUDE
