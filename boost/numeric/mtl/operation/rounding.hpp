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

#ifndef MTL_VEC_ROUNDING_INCLUDE
#define MTL_VEC_ROUNDING_INCLUDE

namespace mtl { namespace vec {


    /// Element-wise ceil of \a v
    template <typename Vector>
    ceil_view<Vector> ceil(const Vector& v)
    {
        return ceil_view<Vector>(v);
    }

    /// Element-wise floor of \a v
    template <typename Vector>
    floor_view<Vector> floor(const Vector& v)
    {
        return floor_view<Vector>(v);
    }

# ifdef MTL_WITH_MATH_ELEVEN    
    /// Element-wise round of \a v
    template <typename Vector>
    round_view<Vector> round(const Vector& v)
    {
        return round_view<Vector>(v);
    }
# endif
    
# ifdef MTL_WITH_MATH_ELEVEN    
    /// Element-wise truncation of \a v
    template <typename Vector>
    trunc_view<Vector> trunc(const Vector& v)
    {
        return trunc_view<Vector>(v);
    }
# endif
    
}} // namespace mtl::vec

#endif // MTL_VEC_ROUNDING_INCLUDE
