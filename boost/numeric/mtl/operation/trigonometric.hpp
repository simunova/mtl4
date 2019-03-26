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

#ifndef MTL_VEC_TRIGONOMETRIC_INCLUDE
#define MTL_VEC_TRIGONOMETRIC_INCLUDE

#include <boost/numeric/mtl/vector/map_view.hpp>

// vector functions

namespace mtl { namespace vec {

    /// Element-wise acos of \a v
    template <typename Vector>
    acos_view<Vector> acos(const Vector& v)
    {
        return acos_view<Vector>(v);
    }

    /// Element-wise acosh of \a v
    template <typename Vector>
    acosh_view<Vector> acosh(const Vector& v)
    {
        return acosh_view<Vector>(v);
    }

    /// Element-wise asin of \a v
    template <typename Vector>
    asin_view<Vector> asin(const Vector& v)
    {
        return asin_view<Vector>(v);
    }

    /// Element-wise asinh of \a v
    template <typename Vector>
    asinh_view<Vector> asinh(const Vector& v)
    {
        return asinh_view<Vector>(v);
    }

    /// Element-wise atan of \a v
    template <typename Vector>
    atan_view<Vector> atan(const Vector& v)
    {
        return atan_view<Vector>(v);
    }

    /// Element-wise atanh of \a v
    template <typename Vector>
    atanh_view<Vector> atanh(const Vector& v)
    {
        return atanh_view<Vector>(v);
    }

    // non-inverse
    
    /// Element-wise cos of \a v
    template <typename Vector>
    cos_view<Vector> cos(const Vector& v)
    {
        return cos_view<Vector>(v);
    }

    /// Element-wise cosh of \a v
    template <typename Vector>
    cosh_view<Vector> cosh(const Vector& v)
    {
        return cosh_view<Vector>(v);
    }

    /// Element-wise sin of \a v
    template <typename Vector>
    sin_view<Vector> sin(const Vector& v)
    {
        return sin_view<Vector>(v);
    }

    /// Element-wise sinh of \a v
    template <typename Vector>
    sinh_view<Vector> sinh(const Vector& v)
    {
        return sinh_view<Vector>(v);
    }

    /// Element-wise tan of \a v
    template <typename Vector>
    tan_view<Vector> tan(const Vector& v)
    {
        return tan_view<Vector>(v);
    }

    /// Element-wise tanh of \a v
    template <typename Vector>
    tanh_view<Vector> tanh(const Vector& v)
    {
        return tanh_view<Vector>(v);
    }


}} // namespace mtl::vec

#endif // MTL_VEC_TRIGONOMETRIC_INCLUDE
