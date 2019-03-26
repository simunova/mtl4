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

#ifndef MTL_SIGNUM_INCLUDE
#define MTL_SIGNUM_INCLUDE

#include <complex>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/operation/real.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl {

namespace sfunctor {

    template <typename Value>
    struct signum
    {
        typedef Value result_type;

        static inline Value apply(const Value& v)
        {
            using math::zero; using math::one;
            return v == zero(v) ? zero(v) : ( v < zero(v) ? -one(v) : one(v) );
        }

        Value operator()(const Value& v) const
        {
            return apply(v);
        }
    };

    template <typename Value>
    struct signum<std::complex<Value> >
    {
        typedef Value result_type;

        static inline Value apply(const std::complex<Value>& v)
        {
            return signum<Value>::apply(mtl::real(v));
        }

        Value operator()(const Value& v) const
        {
            return apply(v);
        }

    };

} // namespace sfunctor

namespace vec {
    
    /// Signum (sign) values of \a v element-wise
    template <typename Vector>
    signum_view<Vector> signum(const Vector& v)
    {
        return signum_view<Vector>(v);
    }
    
} // namespace vec

/// Sign of scalars
/** For complex numbers, the sign of real part is returned; subject to revision. **/
template <typename Value>
inline typename sfunctor::signum<Value>::result_type signum(const Value& v)
{
    vampir_trace<6> tracer;
    return sfunctor::signum<Value>::apply(v);
}


} // namespace mtl

#endif // MTL_SIGNUM_INCLUDE
