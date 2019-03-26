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

#ifndef MTL_FOR_EACH_NONZERO_INCLUDE
#define MTL_FOR_EACH_NONZERO_INCLUDE

#include <utility>

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace vec {

	/// Perform \p f(v[i]) on each non-zero i in constant vector \p v; thus the must keep the result in its state 
	template <typename Vector, typename Functor>
	inline void look_at_each_nonzero(const Vector& v, Functor& f)
	{
	    vampir_trace<2007> tracer;
	    typedef typename mtl::traits::range_generator<tag::const_iter::nz, Vector>::type iterator;
	    for (iterator i= begin<tag::iter::nz>(v), iend= end<tag::iter::nz>(v); i != iend; ++i)
		f(*i);
	}

	/// Perform \p f(v[i], i) on each non-zero i in constant vector \p v; thus the must keep the result in its state 
	template <typename Vector, typename Functor>
	typename mtl::traits::enable_if_vector<Vector>::type // to be called for vectors only 
	inline look_at_each_nonzero_pos(const Vector& v, Functor& f)
	{
	    vampir_trace<2008> tracer;
	    typename mtl::traits::index<Vector>::type           index(v); 
	    typename mtl::traits::const_value<Vector>::type     value(v); 

	    typedef typename traits::range_generator<tag::nz, Vector>::type iterator;
	    for (iterator i= begin<tag::nz>(v), iend= end<tag::nz>(v); i != iend; ++i)
		f(value(*i), index(*i));
	}

    } // namespace vector

    namespace mat {

	/// Perform a potentially mutating \p f(A[i][j]) on each non-zero entry in matrix \p A 
	template <typename Matrix, typename Functor>
	inline void look_at_each_nonzero(const Matrix& A, Functor& f)
	{
	    vampir_trace<3016> tracer;
	    typename mtl::traits::const_value<Matrix>::type     value(A); 

	    typedef typename mtl::traits::range_generator<tag::major, Matrix>::type     cursor_type;
	    typedef typename mtl::traits::range_generator<tag::nz, cursor_type>::type   icursor_type;

	    for (cursor_type cursor = begin<tag::major>(A), cend = end<tag::major>(A); cursor != cend; ++cursor) 
		for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); 
		     icursor != icend; ++icursor)
		    f(value(*icursor));
	}

	/// Perform a potentially mutating \p f(A[i][j], make_pair(i, j)) on each non-zero entry in matrix \p A 
	template <typename Matrix, typename Functor>
	typename mtl::traits::enable_if_matrix<Matrix>::type // to be called for matrices only
	inline look_at_each_nonzero_pos(const Matrix& A, Functor& f)
	{
	    vampir_trace<3017> tracer;
	    typename mtl::traits::row<Matrix>::type             row(A); 
	    typename mtl::traits::col<Matrix>::type             col(A); 
	    typename mtl::traits::const_value<Matrix>::type     value(A); 

	    typedef typename mtl::traits::range_generator<tag::major, Matrix>::type     cursor_type;
	    typedef typename mtl::traits::range_generator<tag::nz, cursor_type>::type   icursor_type;

	    for (cursor_type cursor = begin<tag::major>(A), cend = end<tag::major>(A); cursor != cend; ++cursor) 
		for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); 
		     icursor != icend; ++icursor)
		    f(value(*icursor), std::make_pair(row(*icursor), col(*icursor)));
	}

    } // namespace matrix

} // namespace mtl

#endif // MTL_FOR_EACH_NONZERO_INCLUDE
