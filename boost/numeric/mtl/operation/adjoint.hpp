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

#ifndef MTL_ADJOINT_INCLUDE
#define MTL_ADJOINT_INCLUDE

#include <boost/numeric/mtl/operation/hermitian.hpp>

namespace mtl { namespace mat {

namespace traits {

    template <typename LinOp>
    struct adjoint
    {
	typedef typename mtl::mat::detail::hermitian<LinOp>::result_type  type;
	type operator()(const LinOp& A)
	{
	    return hermitian(A);
	}
    };
}

/// Adjoint linear operator, typically Hermitian transposed
template <typename LinOp>
typename traits::adjoint<LinOp>::type 
inline adjoint(const LinOp& A)
{
    return traits::adjoint<LinOp>()(A);
}

}} // namespace mtl::matrix

#endif // MTL_ADJOINT_INCLUDE
