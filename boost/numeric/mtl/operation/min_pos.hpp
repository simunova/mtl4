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

#ifndef MTL_MIN_POS_INCLUDE
#define MTL_MIN_POS_INCLUDE

#include <utility>

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/pos_type.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/look_at_each_nonzero.hpp>


namespace mtl {

    namespace vec {
	
	template <typename Vector>
	struct min_pos_functor
	{
	    typedef typename Collection<Vector>::value_type       value_type;
	    typedef typename mtl::traits::pos_type<Vector>::type  pos_type;
	    typedef std::pair<value_type, pos_type>               result_type;

	    // initialize with max value and max position
	    min_pos_functor() : value(math::identity(math::min<result_type>(), result_type())) {} 

	    void operator()(const value_type& x, const pos_type& p)
	    {
		if (x < value.first)
		    value= std::make_pair(x, p);
 	    }

	    bool unchanged() const { return value.second == math::identity(math::min<pos_type>(), pos_type()); }

	    result_type  value;
	};
	///Returns position of minimal entry of %vector v
	template <typename Vector>
	typename min_pos_functor<Vector>::pos_type
	inline min_pos(const Vector& v)
	{
	    min_pos_functor<Vector> f;
	    look_at_each_nonzero_pos(v, f);

	    MTL_DEBUG_THROW_IF(f.unchanged(), runtime_error("min_pos cannot be applied on empty container"));
	    return f.value.second;
	}

    } // namespace vector

    namespace mat {

	using mtl::vec::min_pos;
    }


} // namespace mtl

#endif // MTL_MIN_POS_INCLUDE
