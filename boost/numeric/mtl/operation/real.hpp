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

#ifndef MTL_REAL_INCLUDE
#define MTL_REAL_INCLUDE

#include <complex>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl {

namespace sfunctor {

    template <typename Value>
    struct real
    {
	typedef Value result_type;

	static inline Value apply(const Value& v)
	{
	    return v;
	}

	Value operator()(const Value& v) const
	{
	    return v;
	}
    };

    template <typename Value>
    struct real<std::complex<Value> >
    {
	typedef Value result_type;

	static inline result_type apply(const std::complex<Value>& v)
	{
	    return std::real(v);
	}
	result_type operator()(const std::complex<Value>& v) const
	{
	    return std::real(v);
	}
    };
}

/// real part of scalars (including non-complex)
template <typename Value>
typename mtl::traits::enable_if_scalar<Value, typename sfunctor::real<Value>::result_type>::type
inline real(const Value& v)
{	
    vampir_trace<3> tracer;
    return sfunctor::real<Value>::apply(v);
}

namespace vec {

    /// Real part of a vector
    template <typename Vector>
    typename mtl::traits::enable_if_vector<Vector, real_view<Vector> >::type
    inline real(const Vector& v)
    {
	return real_view<Vector>(v);
    }
} 

namespace mat {

    /// Real part of a matrix 
    template <typename Matrix>
    typename mtl::traits::enable_if_matrix<Matrix, real_view<Matrix> >::type
    inline real(const Matrix& v)
    {
	return real_view<Matrix>(v);
    }
} 

} // namespace mtl

#endif // MTL_REAL_INCLUDE
