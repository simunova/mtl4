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

// With contributions from Cornelius Steinhardt

#ifndef MTL_VECTOR_SECULAR_INCLUDE
#define MTL_VECTOR_SECULAR_INCLUDE

#include <cmath>
#include <boost/utility.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/operation/minimal_increase.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace vec {

/// Class for the secular equation( to solve eigenvalue problems)
template <typename Vector>
class secular_f
{
    typedef typename Collection<Vector>::value_type   value_type;
    typedef typename Collection<Vector>::size_type    size_type;

  public:
    /// Constructor needs 2 Vectors z(numerator), d(denominator) and sigma as factor before the sum
    secular_f(const Vector& z, const Vector& d, value_type sigma) 
      : z(z), d(d), sigma(sigma) {}

    /// secular_f equation as function, evaluates the function value
    /** \f$f(x)=1+\sigma * sum_{i=1}^{n}\frac{z_i}{d_i-x} \f$**/
    value_type f(const value_type& lamb)
    {
	value_type fw= 1;
	for(size_type i=0; i<size(z); i++)
	    fw+= sigma*z[i]*z[i]/(d[i]-lamb);
	return fw;
    }

    value_type square(value_type x) const { return x*x; }

    /// gradient of secular_f equation as function, evaluates the gradient function value
    /** \f$gradf(x)=\sigma * sum_{i=1}^{n}\frac{z_i}{(d_i-x)^2} \f$**/
    value_type grad_f(const value_type& lamb)
    {
	value_type gfw= 0.0;
	for(size_type i=0; i<size(z); i++)
	    gfw+= square(z[i] / (d[i] - lamb)); // , std::cout << "gfw = " << gfw << '\n';  //TODO
	return sigma*gfw;
    }
    
    /// Evaluates the roots of secular_f equation =0 with newton algo.
    /** Computes mixed Newton and interval nesting. d must be sorted. **/
    Vector roots()
    {
	assert(size(z) > 1);
	const double tol= 1.0e-6;
	Vector       start(resource(z)), lambda(resource(z));

	for (size_type i= 0; i < size(z); i++) {
	    // Equal poles -> eigenvalue 
	    if (i < size(z) - 1 && d[i] == d[i+1]) { 
		lambda[i]= d[i]; continue; }
	    
	    // Check if root is too close to pole (i.e. d[i]+eps > 0) then take this because we can't reach the root 
	    value_type next= minimal_increase(d[i]), lamb, old;
	    if (f(next) >= value_type(0)){ 
		lambda[i]= next; continue; }
		
	    if (i < size(z) - 1)
		old= lamb= start[i]= (d[i] + d[i+1]) / 2;  //start points between pols
	    else
		old= lamb= start[i]= 1.5 * d[i] - 0.5 * d[i-1];  // last start point plus half the distance to second-last

   	    while (std::abs(f(lamb)) > tol) {
		if (lamb <= d[i])		   
		    start[i]= lamb= (d[i] + start[i]) / 2;  
		else 
		    lamb-= f(lamb) / grad_f(lamb);
		if (old == lamb) break;
		old= lamb;
	    }
	    lambda[i]= lamb;
	} 
	return lambda;
    }

 private:
    Vector     z, d;
    value_type sigma;
};

template <typename Vector, typename Value>
inline Vector secular(const Vector& z, const Vector& d, Value sigma)
{	
	vampir_trace<3030> tracer;
    secular_f<Vector> functor(z, d, sigma);
    return functor.roots();
}

}}// namespace vector


#endif // MTL_VECTOR_SECULAR_INCLUDE

