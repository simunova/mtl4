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

#ifndef MTL_VEC_LOG_INCLUDE
#define MTL_VEC_LOG_INCLUDE

#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl { namespace vec {

    /// Logarithm of \a v element-wise
    template <typename Vector>
    log_view<Vector> log(const Vector& v)
    {
        return log_view<Vector>(v);
    }

# ifdef MTL_WITH_MATH_ELEVEN    
    /// Binary logarithm of \a v element-wise
    template <typename Vector>
    log2_view<Vector> log2(const Vector& v)
    {
        return log2_view<Vector>(v);
    }
# endif

    /// Decimal logarithm of \a v element-wise
    template <typename Vector>
    log10_view<Vector> log10(const Vector& v)
    {
        return log10_view<Vector>(v);
    }


}} // namespace mtl::vec

#endif // MTL_VEC_LOG_INCLUDE
