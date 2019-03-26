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

#ifndef MTL_VEC_SQRT_INCLUDE
#define MTL_VEC_SQRT_INCLUDE

#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl { namespace vec {

    /// Square root of \a v element-wise
    template <typename Vector>
    sqrt_view<Vector> sqrt(const Vector& v)
    {
        return sqrt_view<Vector>(v);
    }

    /// Inverse square root of \a v element-wise
    template <typename Vector>
    rsqrt_view<Vector> rsqrt(const Vector& v)
    {
        return rsqrt_view<Vector>(v);
    }


}} // namespace mtl::vec

#endif // MTL_VEC_SQRT_INCLUDE
