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

#ifndef MTL_MATRIX_MAX_ABS_POS_INCLUDE
#define MTL_MATRIX_MAX_ABS_POS_INCLUDE

#include <utility>
#include <cmath>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/concept/magnitude.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { 

namespace mat {
///Returns pair (row, col) from absolut maximal entry of %matrix A
    template <typename Matrix>
    typename mtl::traits::enable_if_matrix<Matrix, std::pair<typename Collection<Matrix>::size_type, typename Collection<Matrix>::size_type> >::type
    inline max_abs_pos(const Matrix& A)
    {
	vampir_trace<3024> tracer;
	namespace traits = mtl::traits;
	using std::abs;
	typedef typename Collection<Matrix>::value_type   value_type;
	typedef typename Collection<Matrix>::size_type    size_type;

	typename RealMagnitude<value_type>::type max(abs(A[0][0]));
	size_type r= 0, c= 0;

	typename traits::row<Matrix>::type             row(A); 
	typename traits::col<Matrix>::type             col(A); 
	typename traits::const_value<Matrix>::type     value(A); 
	typedef typename traits::range_generator<tag::major, Matrix>::type  cursor_type;
	
	for (cursor_type cursor = begin<tag::major>(A), cend = end<tag::major>(A); cursor != cend; ++cursor) {
	    typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	    for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) 
		if (abs(value(*icursor)) > max) {
		    max= abs(value(*icursor));
		    r= row(*icursor);
		    c= col(*icursor);
		}
	}
	
	return std::make_pair(r, c);
    }

} // namespace matrix

namespace vec {
///Returns position from absolut maximal entry of %vector v
    template <typename Vector>
    typename mtl::traits::enable_if_vector<Vector, typename Collection<Vector>::size_type>::type
    inline max_abs_pos(const Vector& v)
    {
	vampir_trace<2011> tracer;
	using std::abs;
	typedef typename Collection<Vector>::size_type    size_type;
	typedef typename Collection<Vector>::value_type   value_type;

	size_type i= 0;
	size_type max_col= size(v);
	value_type max(abs(v[0]));
	
	for(size_type j= 1; j < max_col; j++)
	    if(abs(v[j]) > max) {
		max = abs(v[j]);
		i= j;
	    }
	return i;
    }

} // namespace vector

} // namespace mtl

#endif // MTL_MATRIX_MAX_ABS_POS_INCLUDE

