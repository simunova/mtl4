// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling, Cornelius Steinhardt and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef ITL_ARMIJO_INCLUDE
#define ITL_ARMIJO_INCLUDE

namespace itl {

/// Step size control by Armijo
/** 
 **/
template <typename Value= double>
class armijo
{
  public:
    typedef Value   value_type;

    // Defaults from Prof. Fischer's lecture
    armijo(Value delta= 0.5, Value gamma= 0.5, Value beta1= 0.25, Value beta2= 0.5)
      : delta(delta), gamma(gamma), beta1(beta1), beta((beta1 + beta2) / 2.0) {}

    /// 
    template <typename Vector, typename F, typename Grad>
    typename mtl::Collection<Vector>::value_type 
    operator() (const Vector& x, const Vector& d, F f, Grad grad_f) const
    {
	// Star's step size
	typename mtl::Collection<Vector>::value_type alpha= -gamma * dot(grad_f(x), d) / dot(d, d);
	Vector     x_k(x + alpha * d);

	while (f(x_k) > f(x) + (beta1 * alpha) * dot(grad_f(x), d)) {	
	    alpha*= beta;
	    x_k= x+ alpha * d;
	}
	return alpha;
    }
  private:
    Value delta, gamma, beta1, beta;
};


} // namespace itl

#endif // ITL_ARMIJO_INCLUDE
