// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_STATIC_VECTOR_INCLUDE
#define MTL_STATIC_VECTOR_INCLUDE

#if defined(MTL_WITH_STATICASSERT) && defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)

#include <cstddef>

namespace mtl {

    namespace impl {

	// Need helper, cannot specialize class inside class
	template <typename Vector, std::size_t Pos>
	struct static_vector_get
	{
	    static_assert(Pos < Vector::length, "Position is out of range.");
	    typedef typename static_vector_get<typename Vector::tail, Pos - 1>::type type;
	};

	template <typename Vector>
	struct static_vector_get<Vector, 0>
	{
	    static_assert(Vector::length >= 1, "Position is out of range.");
	    typedef Vector   type;
	};

    }


    template <typename ValueType, ValueType FirstValue, ValueType ...Values>
    struct static_vector
    {
	typedef ValueType                                        value_type;
	typedef static_vector<ValueType, FirstValue, Values...>  self;

	static const std::size_t                                 length= sizeof...(Values) + 1;
	static const ValueType                                   value= FirstValue;

	typedef static_vector<ValueType, Values...>              tail;

	template <std::size_t Pos>
	using get= typename impl::static_vector_get<self, Pos>::type;
    };

    template <typename ValueType, ValueType FirstValue>
    struct static_vector<ValueType, FirstValue>
    {
	typedef ValueType                                        value_type;
	typedef static_vector<ValueType, FirstValue>             self;

	static const std::size_t                                 length= 1;
	static const ValueType                                   value= FirstValue;

	template <std::size_t Pos>
	using get= typename impl::static_vector_get<self, Pos>::type;
    };

} // namespace mtl

#endif // required C++11 features

#endif // MTL_STATIC_VECTOR_INCLUDE
