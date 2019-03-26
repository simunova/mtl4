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

#ifndef MTL_CONJ_INCLUDE
#define MTL_CONJ_INCLUDE

#include <complex>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/algebraic_category.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>
#include <boost/numeric/mtl/utility/view_code.hpp>
#include <boost/numeric/mtl/utility/viewed_collection.hpp>
#include <boost/numeric/mtl/utility/compose_view.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/matrix/view_ref.hpp>
#include <boost/numeric/mtl/matrix/map_view.hpp>
#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl {

namespace sfunctor {

    template <typename Value, typename AlgebraicCategory>
    struct conj_aux
    {
	typedef Value result_type;

	static inline result_type apply(const Value& v)
	{
	    return v;
	}

	result_type operator() (const Value& v) const
	{
	    return v;
	}
    };


    template <typename Value, typename AlgebraicCategory>
    struct conj_aux<std::complex<Value>, AlgebraicCategory>
    {
	typedef std::complex<Value> result_type;

	static inline result_type apply(const std::complex<Value>& v)
	{
	    return std::conj(v);
	}

	result_type operator() (const std::complex<Value>& v) const
	{
	    return std::conj(v);
	}
    };

    // Only declarations here, definitions in mat::map_view (map_view)
    template <typename Matrix>
    struct conj_aux<Matrix, tag::matrix>;

    template <typename Vector>
    struct conj_aux<Vector, tag::vector>;

    // Short cut for result type
    template <typename Value>
    struct conj
	: public conj_aux<Value, typename mtl::traits::algebraic_category<Value>::type>
    {};

} // namespace sfunctor
    
    namespace vec {

	/// Conjugate of an vector
	template <typename Vector>
	typename mtl::traits::enable_if_vector<Vector, conj_view<Vector> >::type
	inline conj(const Vector& v)
	{
	    return conj_view<Vector>(v);
	}
    } 

    namespace mat {

	namespace detail {

	    template <typename Matrix>
	    struct conj_trait
	    {
		static const unsigned code= mtl::traits::view_toggle_conj<mtl::traits::view_code<Matrix> >::value;
		typedef typename mtl::traits::compose_view<code, typename mtl::traits::viewed_collection<Matrix>::type>::type type;
		
		static inline type apply(const Matrix& A)
		{
		    return type(view_ref(A));
		}
	    };

	}

	/// Conjugate of a matrix
	template <typename Matrix>
	typename mtl::traits::lazy_enable_if_matrix<Matrix, detail::conj_trait<Matrix> >::type
	inline conj(const Matrix& A)
	{
	    return detail::conj_trait<Matrix>::apply(A);
	}
    } 

    namespace scalar {

	// Only scalar values remain here
	template <typename Value>
	typename mtl::traits::enable_if_scalar<
	    Value
	  , typename sfunctor::conj<Value>::result_type
	>::type
	inline conj(const Value& v)
	{
	    return mtl::sfunctor::conj<Value>::apply(v);
	}

	float inline conj(float v) { return v; }
	double inline conj(double v) { return v; }
	long double inline conj(long double v) { return v; }
    }

    /// Conjugate of vector, matrix, or scalar
    using vec::conj;
    using mat::conj; 
    using scalar::conj;
    // using std::conj;

} // namespace mtl

#endif // MTL_CONJ_INCLUDE
