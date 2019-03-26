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

#ifndef ITL_FSM_INCLUDE
#define ITL_FSM_INCLUDE

#include <boost/numeric/mtl/operation/two_norm.hpp>

namespace itl {

/// Folded spectrum method
/** Computed and named as in http://en.wikipedia.org/wiki/Folded_spectrum_method **/
template < typename LinearOperator, typename VectorSpace, typename EigenValue, 
	   typename Damping, typename Iteration >
int fsm(const LinearOperator& H, VectorSpace& phi, EigenValue eps, Damping alpha, Iteration& iter)
{
    VectorSpace v1(H * phi - eps * phi);
    for (; !iter.finished(v1); ++iter) {
	VectorSpace v2(H * v1 - eps * v1);
	phi-= alpha * v2;
	phi/= two_norm(phi);
	v1= H * phi - eps * phi;
    }
    return iter;
}


} // namespace itl

#endif // ITL_FSM_INCLUDE
