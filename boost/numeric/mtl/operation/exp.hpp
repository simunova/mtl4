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

#ifndef MTL_VEC_EXP_INCLUDE
#define MTL_VEC_EXP_INCLUDE

#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>

namespace mtl { 

    namespace vec {

	/// Exponential of \a v element-wise
	template <typename Vector>
	exp_view<Vector> exp(const Vector& v)
	{
	    return exp_view<Vector>(v);
	}

#      ifdef MTL_WITH_MATH_ELEVEN    
	/// Binary exponential of \a v element-wise
	template <typename Vector>
	exp2_view<Vector> exp2(const Vector& v)
	{
	    return exp2_view<Vector>(v);
	}
#      endif

	/// Decimal exponential of \a v element-wise
	template <typename Vector>
	exp10_view<Vector> exp10(const Vector& v)
	{
	    return exp10_view<Vector>(v);
	}

    } // namespace vec

    namespace mat {

	/// Exponential of \a A element-wise
	template <typename Matrix>
	exp_view<Matrix> exp(const Matrix& A)
	{
	    return exp_view<Matrix>(A);
	}

    } // namespace mat


} // namespace mtl

#endif // MTL_VEC_EXP_INCLUDE
