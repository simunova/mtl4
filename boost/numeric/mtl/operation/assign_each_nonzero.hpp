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

#ifndef MTL_ASSIGN_EACH_NONZERO_INCLUDE
#define MTL_ASSIGN_EACH_NONZERO_INCLUDE

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

    namespace vec {

	/// Assign the result of \p f(i) to each non-zero i in non-constant vector \p v
	template <typename Vector, typename Functor>
	inline void assign_each_nonzero(Vector& v, const Functor& f)
	{
	    vampir_trace<3008> tracer;
	    typedef typename traits::range_generator<tag::iter::nz, Vector>::type iterator;
	    for (iterator i= begin<tag::iter::nz>(v), iend= end<tag::iter::nz>(v); i != iend; ++i)
		*i= f(*i);
	}

    } // namespace vector

    namespace mat {

	/// Assign the result of \p f(i) to each non-zero i in non-constant matrix \p A 
	template <typename Matrix, typename Functor>
	inline void assign_each_nonzero(Matrix& m, const Functor& f)
	{
	    vampir_trace<3008> tracer;
	    typename mtl::traits::value<Matrix>::type     value(m); 

	    typedef typename mtl::traits::range_generator<tag::major, Matrix>::type     cursor_type;
	    typedef typename mtl::traits::range_generator<tag::nz, cursor_type>::type   icursor_type;

	    for (cursor_type cursor = begin<tag::major>(m), cend = end<tag::major>(m); cursor != cend; ++cursor) 
		for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); 
		     icursor != icend; ++icursor)
			value(*icursor, f(value(*icursor))); 
	}

    } // namespace matrix

} // namespace mtl

#endif // MTL_ASSIGN_EACH_NONZERO_INCLUDE
