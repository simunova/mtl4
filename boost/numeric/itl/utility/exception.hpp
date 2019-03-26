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

#ifndef ITL_EXCEPTION_INCLUDE
#define ITL_EXCEPTION_INCLUDE

#include <boost/numeric/mtl/utility/exception.hpp>

namespace itl {

/// Exception for iterative solvers that exhausted the search space, i.e. search direction(s) parallel to already visited Krylov subspace
/** Either your matrix is too badly conditioned or your termination criterion is too strict for the used arithmetic (maybe use Gnu Multi-precision library). **/
struct search_space_exhaustion
  : mtl::runtime_error
{
    /// Error can be specified more precisely in constructor if desired
    explicit search_space_exhaustion(const char *s= "Iterative solvers that exhausted the search space, i.e. search direction(s) parallel to already visited Krylov subspace")
      : mtl::runtime_error(s) {}
};

/// Vectors unexpectedly become orthogonal, special case of search_space_exhaustion.
struct unexpected_orthogonality : search_space_exhaustion {};

} // namespace itl

#endif // ITL_EXCEPTION_INCLUDE
