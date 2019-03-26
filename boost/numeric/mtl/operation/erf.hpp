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

#ifndef MTL_VEC_ERF_INCLUDE
#define MTL_VEC_ERF_INCLUDE

#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl { namespace vec {

# ifdef MTL_WITH_MATH_ELEVEN    
    /// Error function of \a v element-wise
    template <typename Vector>
    erf_view<Vector> erf(const Vector& v)
    {
        return erf_view<Vector>(v);
    }

    /// Inverse error function of \a v element-wise
    template <typename Vector>
    erfc_view<Vector> erfc(const Vector& v)
    {
        return erfc_view<Vector>(v);
    }
# endif

}} // namespace mtl::vec

#endif // MTL_VEC_ERF_INCLUDE
