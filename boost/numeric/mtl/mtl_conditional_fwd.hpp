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

#ifndef MTL_MTL_CONDITIONAL_FWD_INCLUDE
#define MTL_MTL_CONDITIONAL_FWD_INCLUDE

// Forward declarations that need meta-programming (enable_if and alike)

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_composable_vector.hpp>
#include <boost/numeric/mtl/utility/is_distributed.hpp>
#include <boost/numeric/mtl/utility/iset.hpp>
#include <boost/numeric/mtl/utility/is_lazy.hpp>
#include <boost/numeric/mtl/utility/is_multi_vector_expr.hpp>
#include <boost/numeric/mtl/utility/is_row_major.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>
#include <boost/numeric/mtl/utility/is_vector_reduction.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>


namespace mtl {

    namespace vec {

	template <typename Vector>
	typename mtl::traits::enable_if_vector<Vector, conj_view<Vector> >::type
	inline conj(const Vector& v);
    }

} // namespace mtl

#endif // MTL_MTL_CONDITIONAL_FWD_INCLUDE
