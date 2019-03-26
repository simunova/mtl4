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

#ifndef ITL_PSB_INCLUDE
#define ITL_PSB_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/itl/utility/exception.hpp>

namespace itl {

/// Update of Hessian matrix for e.g. Quasi-Newton by Powell's symmetric Broyden formula
struct psb
{
    /// \f$ H_{k+1}=B_{k+1}^{-1}=H_k+\frac{(y_k-H_k\cdot s_k)s_k^T+s_k(y_k-H_k\cdot s_k)^T}{s_k^T\cdot s_k}-\frac{(y_k-H_k\cdot s_k)^T\cdot s_k}{(s_k^ts_k)^2}s_k\cdot s_k^T \f$
    template <typename Matrix, typename Vector>
    void operator() (Matrix& H, const Vector& y, const Vector& s)
    {
	typedef typename mtl::Collection<Vector>::value_type value_type;
	assert(num_rows(H) == num_cols(H));
	Vector     a(s - H * y);
	value_type gamma= 1 / dot (y, y);
        MTL_THROW_IF(gamma == 0.0, unexpected_orthogonality());
    
        H+= gamma * a * trans(y) + gamma * y * trans(a) - dot(a, y) * gamma * gamma * y * trans(y);
   }
};



} // namespace itl

#endif // ITL_PSB_INCLUDE

