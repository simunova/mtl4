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

#ifndef GLAS_GLAS_TAG_INCLUDE
#define GLAS_GLAS_TAG_INCLUDE

#undef major

namespace glas { namespace tag {

// To iterate only over non-zero elements
struct nz {};

// To iterate over all elements
struct all {};

// To iterate over rows
// Generated cursors must provide range generators
struct row {};

// To iterate over cols
// Generated cursors must provide range generators
struct col {};

// To iterate over the major dimension of matrices (like MTL 2)
struct major {};

// Same with minor
struct minor {};

}} // namespace glas::tag

#endif // GLAS_GLAS_TAG_INCLUDE
