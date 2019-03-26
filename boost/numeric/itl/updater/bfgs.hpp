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

#ifndef ITL_BFGS_INCLUDE
#define ITL_BFGS_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/itl/utility/exception.hpp>

namespace itl {

/// Update of Hessian matrix for e.g. Quasi-Newton by Broyden, Fletcher, Goldfarb, and Shanno
struct bfgs
{
    /// \f$ H_{k+1}=B_{k+1}^{-1}=(I-\frac{y_k\cdot s_k^T}{y_k^T\cdot s_k})^T\cdot H_k \cdot (I-\frac{y_k\cdot s_k^T}{y_k^T\cdot s_k}) + \frac{s_k\cdot s_k^T}{y_k^T\cdot s_k}\f$ 
    template <typename Matrix, typename Vector>
    void operator() (Matrix& H, const Vector& y, const Vector& s)
    {
	typedef typename mtl::Collection<Vector>::value_type value_type;
	assert(num_rows(H) == num_cols(H));

	value_type gamma= 1 / dot(y,s);
	MTL_THROW_IF(gamma == 0.0, unexpected_orthogonality());
	Matrix     A(math::one(H) - gamma * s * trans(y)),
	           H2(A * H * trans(A) + gamma * s * trans(s));
	swap(H2, H); // faster than H= H2
    }
}; 



} // namespace itl

#endif // ITL_BFGS_INCLUDE
