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

#ifndef MTL_IMAG_INCLUDE
#define MTL_IMAG_INCLUDE

#include <complex>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>

namespace mtl {

namespace sfunctor {

    template <typename Value>
    struct imag
    {
	typedef Value result_type;

	static inline result_type apply(const Value& v)
	{
	    using math::zero;
	    return zero(v);
	}
	result_type operator()(const Value& v) const
	{
	    using math::zero;
	    return zero(v);
	}
    };

    template <typename Value>
    struct imag<std::complex<Value> >
    {
	typedef Value result_type;

	static inline result_type apply(const std::complex<Value>& v)
	{
	    return std::imag(v);
	}
	result_type operator()(const std::complex<Value>& v) const
	{
	    return std::imag(v);
	}
    };

}

/// Imaginary part of scalars (including non-complex)
template <typename Value>
typename mtl::traits::enable_if_scalar<Value, typename sfunctor::imag<Value>::result_type>::type
inline imag(const Value& v)
{
    return sfunctor::imag<Value>::apply(v);
}

namespace vec {

    /// Imaginary part of an vector
    template <typename Vector>
    typename mtl::traits::enable_if_vector<Vector, imag_view<Vector> >::type
    inline imag(const Vector& v)
    {
	return imag_view<Vector>(v);
    }
} 

namespace mat {

    /// Imaginary part of an vector
    template <typename Matrix>
    typename mtl::traits::enable_if_matrix<Matrix, imag_view<Matrix> >::type
    inline imag(const Matrix& v)
    {
	return imag_view<Matrix>(v);
    }
} 


} // namespace mtl

#endif // MTL_IMAG_INCLUDE
