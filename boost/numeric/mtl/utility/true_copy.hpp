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

#ifndef MTL_TRAITS_TRUE_COPY_INCLUDE
#define MTL_TRAITS_TRUE_COPY_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>

namespace mtl { namespace traits {

/// Type trait to force copy
/** Some libs define types to force shallow copy. This causes stale references in expression templates.
    To counter-act we substitute the types, e.g. mp_real_tmp with mp_real, to force a copy. **/
template <typename T>
struct true_copy
{
    typedef T type;
};



}} // namespace mtl::traits

#endif // MTL_TRAITS_TRUE_COPY_INCLUDE
