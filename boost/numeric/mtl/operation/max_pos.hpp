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

#ifndef MTL_MATRIX_MAX_POS
#define MTL_MATRIX_MAX_POS

#include <cmath>
#include <utility>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/pos_type.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/look_at_each_nonzero.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

#if 0
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#endif



namespace mtl { 

    namespace vec {
	
	template <typename Vector>
	struct max_pos_functor
	{
	    typedef typename Collection<Vector>::value_type       value_type;
	    typedef typename mtl::traits::pos_type<Vector>::type  pos_type;
	    typedef std::pair<value_type, pos_type>               result_type;

	    // initialize with max value and max position
	    max_pos_functor() 
	    {
		value.first= math::identity(math::max<value_type>(), value_type()); // minimal value for comparison
		value.second= math::identity(math::min<pos_type>(), pos_type());    // maximal position to check if changed
	    }

	    void operator()(const value_type& x, const pos_type& p)
	    {
		if (x > value.first)
		    value= std::make_pair(x, p);
	    }

	    bool unchanged() const { return value.second == math::identity(math::min<pos_type>(), pos_type()); }

	    result_type  value;
	};
	///Returns position of maximal entry of %vector v
	template <typename Vector>
	typename max_pos_functor<Vector>::pos_type
	inline max_pos(const Vector& v)
	{
	    vampir_trace<2013> tracer;
	    max_pos_functor<Vector> f;
	    look_at_each_nonzero_pos(v, f);

	    MTL_DEBUG_THROW_IF(f.unchanged(), runtime_error("max_pos cannot be applied on empty container"));
	    return f.value.second;
	}

    } // namespace vector

    namespace mat {
	///Returns pair (row, col) from maximal entry of %matrix A
	using mtl::vec::max_pos;
    }


#if 0
namespace matrix {
///Returns pair (row, col) from maximal entry of %matrix A
    template <typename Matrix>
    std::pair<typename Collection<Matrix>::size_type, typename Collection<Matrix>::size_type>
    inline max_pos(const Matrix& A)
    {
	namespace traits = mtl::traits;
	typedef typename Collection<Matrix>::value_type   value_type;
	typedef typename Collection<Matrix>::size_type    size_type;

	value_type max(A[0][0]);
	size_type r= 0, c= 0;

	typename traits::row<Matrix>::type             row(A); 
	typename traits::col<Matrix>::type             col(A); 
	typename traits::const_value<Matrix>::type     value(A); 
	typedef typename traits::range_generator<tag::major, Matrix>::type  cursor_type;
	
	for (cursor_type cursor = begin<tag::major>(A), cend = end<tag::major>(A); cursor != cend; ++cursor) {
	    typedef typename traits::range_generator<tag::nz, cursor_type>::type icursor_type;
	    for (icursor_type icursor = begin<tag::nz>(cursor), icend = end<tag::nz>(cursor); icursor != icend; ++icursor) 
		if (value(*icursor) > max) {
		    max= value(*icursor);
		    r= row(*icursor);
		    c= col(*icursor);
		}
	}
	
	return std::make_pair(r, c);
}

} // namespace matrix

namespace vector {
///Returns position from maximal entry of %vector v
    template <typename Vector>
    typename Collection<Vector>::size_type
    inline max_pos(const Vector& v)
    {
	typedef typename Collection<Vector>::size_type    size_type;
	typedef typename Collection<Vector>::value_type   value_type;

	size_type i= 0;
	
	size_type max_col= size(v);
	value_type max(v[0]);
	
	for(size_type j= 1;i < max_col; j++)
	    if(v[j] > max) {
		max = v[j];
		i= j;
	    }
	return i;
    }

} // namespace vector

#endif


} // namespace mtl

#endif // MTL_MATRIX_MAX_POS

