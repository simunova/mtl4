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

#ifndef MTL_VECTOR_UNIT_VECTOR_INCLUDE
#define MTL_VECTOR_UNIT_VECTOR_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>

namespace mtl { 

namespace traits {

    /// Type of unit_vector; will be changed later to a proxy for the sake of efficiency
    template  <typename Value= double>
    struct unit_vector
    {
	typedef mtl::vec::dense_vector<Value, mtl::vec::parameters<> >  type;
    };
}

namespace vec {

    /// Return k-th unit vector of size n
    /** The result is a dense column vector. In the future this will
	be replaced by a proxy for the sake of efficiency.
	If you use unit_vector in an expression you will not encounter
	this change. If you define a variable you should use
	traits::unit_vector, e.g.:
	\code
	typename mtl::traits::unit_vector<float>::type  e_k(mtl::unit_vector(k, n));
	\endcode
    **/
    template <typename Value>
    typename traits::unit_vector<Value>::type
    inline unit_vector(std::size_t k, std::size_t n)
    {
	using math::zero; using math::one;
	dense_vector<Value> v(n, zero(Value()));
	v[k]= one(Value());
	return v;
    }

    /// Unit vector of type double
    traits::unit_vector<double>::type
    inline unit_vector(std::size_t k, std::size_t n)
    {
	return unit_vector<double>(k, n);
    }

} // namespace mtl::vector

} // namespace mtl

#endif // MTL_VECTOR_UNIT_VECTOR_INCLUDE
