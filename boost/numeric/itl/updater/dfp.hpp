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

#ifndef ITL_DFP_INCLUDE
#define ITL_DFP_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/itl/utility/exception.hpp>

namespace itl {

/// Update of Hessian matrix for e.g. Quasi-Newton by Davidon, Fletcher and Powell formula
struct dfp
{
    /// \f$ H_{k+1}=B_{k+1}^{-1}=H_k+\frac{s_k\cdot s_k^T}{y_k^T\cdot s_k}- \frac{H_k\cdot y_k\cdot y_k^T\cdot H_k^T}{y_k^TH_k\cdot y_k}\f$
    template <typename Matrix, typename Vector>
    void operator() (Matrix& H, const Vector& y, const Vector& s)
    {
	typedef typename mtl::Collection<Vector>::value_type value_type;
	assert(num_rows(H) == num_cols(H));

	Vector h(H*y);
	value_type gamma= 1 / dot(y,s), alpha= 1 / dot(y,h);
	MTL_THROW_IF(gamma == 0.0, unexpected_orthogonality());
	MTL_THROW_IF(alpha == 0.0, unexpected_orthogonality());
	Matrix     A(alpha * y * trans(y)),
	           H2(H - H * A * H  + gamma * s * trans(s));
	swap(H2, H); // faster than H= H2
   }
}; 



} // namespace itl

#endif // ITL_DFP_INCLUDE
