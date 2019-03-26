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

#ifndef ITL_DEFAULT_ITERATION_INCLUDE
#define ITL_DEFAULT_ITERATION_INCLUDE

#include <boost/numeric/itl/iteration/cyclic_iteration.hpp>

namespace itl {

/// Default cyclic iteration
/** max. 1000 iterations, relative tolerance 1e-10, absolute tolerance 0, set to quite. **/
template <class Real>
inline cyclic_iteration<Real> default_iteration()
{
    cyclic_iteration<Real> i(Real(0), 100, Real(1e-10), Real(0));
    i.set_quite(true);
    return i;
}


} // namespace itl

#endif // ITL_DEFAULT_ITERATION_INCLUDE
