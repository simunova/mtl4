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

#ifndef MTL_MATRIX_HOUSEHOLDER_INCLUDE
#define MTL_MATRIX_HOUSEHOLDER_INCLUDE

#include <cmath>
#include <cassert>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl { namespace vec {


/// Computes Householder vector v and scalar b for vector \p y 
/** such that identity_matrix(size(y))-b*v*v' projects the vector y 
    to a positive multiple of the first unit vector. **/
template <typename Vector>
std::pair<typename mtl::vec::dense_vector<typename Collection<Vector>::value_type, parameters<> >, typename Collection<Vector>::value_type>
inline householder(Vector& y)
{
    vampir_trace<2004> tracer;
    assert(size(y) > 0);
    typedef typename  Collection<Vector>::value_type   value_type;
    // typedef typename  Collection<Vector>::size_type    size_type;
    const value_type  zero= math::zero(y[0]), one= math::one(y[0]);

    Vector            v(y);
    v[0]= one;
    irange            tail(1, imax); 
    value_type        s( dot(v[tail], v[tail]) ), b, v0;

    //evaluation of v and b
    if (s == zero)
        b= zero;
    else {
	value_type mu= sqrt(y[0] * y[0] + s);
	v0= v[0]= y[0] <= zero ? y[0] - mu : -s / (y[0] + mu); // complex < zero????
	b= 2 * v0 * v0 / (s + v0 * v0);                        // 2* complex???????
	v/= v0;                                                // normalization of the first entry
    }
  
    return std::make_pair(v,b);
}

/// Computes Householder vector for vector \p y 
/** More stable Householder transformation, also for non-square matrices. **/
template <typename Vector>
typename mtl::vec::dense_vector<typename Collection<Vector>::value_type, parameters<> >
inline householder_s(Vector& y)
{
    vampir_trace<2005> tracer;
    typedef typename  Collection<Vector>::value_type   value_type;
    const value_type  zero= math::zero(y[0]);

    Vector            u(y);
    value_type        nu(sqrt( dot(u, u) )), s;

    if (nu != zero){
	if(u[0] < 0){
		s= -1;
	} else {
		s= 1; 
	}
	u[0]= u[0] + s * nu;
	u/= sqrt( dot(u, u) );
    }

    return u;
}


}} // namespace mtl::matrix

#endif // MTL_MATRIX_HOUSEHOLDER_INCLUDE

