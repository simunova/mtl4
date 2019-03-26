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

#ifndef MTL_MATRIX_SOLVE_INCLUDE
#define MTL_MATRIX_SOLVE_INCLUDE

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/lu.hpp>
#include <boost/numeric/mtl/interface/umfpack_solve.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {

    namespace detail {

	template <typename Matrix, typename Vector>
	Vector inline solve(const Matrix& A, const Vector& b, tag::dense)
	{
		vampir_trace<3034> tracer;
	    return lu_solve(A, b);
	}

#     ifdef MTL_HAS_UMFPACK
	template <typename Value, typename Parameters, typename Vector>
	Vector inline solve(const Matrix& A, const Vector& b, tag::compressed2D)
	{
		vampir_trace<3035> tracer;
	    Vector x(num_cols(A));
	    umfpack_solve(A, x, b);
	    return x;
	}
#     endif
    }


template <typename Matrix, typename Vector>
Vector inline solve(const Matrix& A, const Vector& b)
{
    return detail::solve(A, b, typename category<Coll>::type());
}

}} // namespace mtl::matrix

#endif // MTL_MATRIX_SOLVE_INCLUDE
