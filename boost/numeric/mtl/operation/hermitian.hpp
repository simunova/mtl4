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

#ifndef MTL_HERMITIAN_INCLUDE
#define MTL_HERMITIAN_INCLUDE

#include <boost/numeric/mtl/matrix/hermitian_view.hpp>
#include <boost/numeric/mtl/matrix/view_ref.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/view_code.hpp>
#include <boost/numeric/mtl/utility/viewed_collection.hpp>
#include <boost/numeric/mtl/utility/compose_view.hpp>

namespace mtl { 

    // vector version to be done

    namespace mat {

	namespace detail {

	    template <typename Matrix>
	    struct hermitian
	    {
		static const unsigned code= mtl::traits::view_toggle_hermitian<mtl::traits::view_code<Matrix> >::value;
		typedef typename mtl::traits::compose_view<code, typename mtl::traits::viewed_collection<Matrix>::type>::type result_type;
	
		static inline result_type apply(const Matrix& A)
		{
		    return result_type(view_ref(A));
		}
	    };

	} // namespace detail
	
	/// Return hermitian of matrix A
	template <typename Matrix>
	typename mtl::traits::enable_if_matrix<Matrix, typename detail::hermitian<Matrix>::result_type >::type
	inline hermitian(const Matrix& A)
	{
	    return detail::hermitian<Matrix>::apply(A);
	}
    }

    using mat::hermitian;

} // namespace mtl

#endif // MTL_HERMITIAN_INCLUDE
