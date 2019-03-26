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

#ifndef MTL_VECTOR_TYPE_GENERATOR_INCLUDE
#define MTL_VECTOR_TYPE_GENERATOR_INCLUDE

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#include <cstddef>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/type_parameter.hpp>
#include <boost/numeric/mtl/vector/dimension.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>

namespace mtl {


    namespace type_para {

	template <typename DimPara>
	struct set_vector_dimensions
	{
	    static_assert(! boost::is_same<DimPara, DimPara>::value, "Unsupported argument for vector dimension.");
	};

	template <>
	struct set_vector_dimensions<none>
	{
	    typedef mtl::vec::non_fixed::dimension type;
	};

# ifndef _MSC_VER // creates problems on VS
	template <std::size_t ...Values>
	struct set_vector_dimensions<dim<Values...> >
	{
	    static_assert(sizeof...(Values) == 1, "dim<size> must have exactly 1 argument for vectors!");	    
	};
# endif

	template <std::size_t Size>
	struct set_vector_dimensions<dim<Size> >
	{
	    typedef mtl::vec::fixed::dimension<Size> type;
	};

	
	template <typename TypePara>
	struct vector_parameter_generator
	{
	    typedef typename boost::mpl::at<TypePara, orientation>::type            ori1;
	    typedef typename boost::mpl::if_<boost::is_same<ori1, none>,
					     col_major,
					     ori1>::type                            ori2;

	    typedef typename set_vector_dimensions<typename boost::mpl::at<TypePara, dimensionality>::type>::type dim_type;
	    typedef typename boost::mpl::at<TypePara, size_type>::type                                            as_size;
	    typedef mtl::vec::parameters<
		ori2,
		dim_type,
		boost::is_same<typename boost::mpl::at<TypePara, location>::type, on_stack>::value,
		typename as_size::type
	    > type;
	};


	template <typename Value, typename TypePara>
	struct vector_type_generator
	{
	    typedef typename vector_parameter_generator<TypePara>::type  vector_parameters;
	    typedef typename boost::mpl::if_<
		boost::is_same<typename boost::mpl::at<TypePara, density>::type, sparse>,
		mtl::vec::sparse_vector<Value, vector_parameters>,
		mtl::vec::dense_vector<Value, vector_parameters>
		>::type                                                  type;
	};

    } // namespace type_para

# ifdef MTL_WITH_TEMPLATE_ALIAS
    template <typename Value, typename ...Parameters>
    using vector= typename type_para::vector_type_generator<Value, typename set_parameters<Parameters...>::type>::type;
# endif

} // namespace mtl

#endif // defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#endif // MTL_VECTOR_TYPE_GENERATOR_INCLUDE
