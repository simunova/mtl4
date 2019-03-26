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

#ifndef MTL_TYPE_PARAMETER_COMMON_INCLUDE
#define MTL_TYPE_PARAMETER_COMMON_INCLUDE

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#include <utility>
#include <cstddef>

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/void.hpp>
#include <boost/mpl/if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>


namespace mtl {

    namespace type_para {

	// Helper type for still unset parameter
	struct none {};

	// Parameter kinds (extensible)
	struct density {};
	struct layout {};
	struct symmetry {};
	struct orientation {};
	struct dimensionality {};
	struct size_type {};
	struct location {};
	struct masking {};
	
	// Initial map for common parameter kinds (extensible)
	typedef boost::mpl::map<
	    boost::mpl::pair<density, none>, 
	    boost::mpl::pair<layout, none>, 
	    boost::mpl::pair<symmetry, none>, 
	    boost::mpl::pair<orientation, none>, 
	    boost::mpl::pair<dimensionality, none>, 
	    boost::mpl::pair<size_type, none>, 
	    boost::mpl::pair<location, none>,
	    boost::mpl::pair<masking, none>
	>                            init_common;

	// The kinds of common parameters except templated parameters
	typedef boost::mpl::map<
	    boost::mpl::pair<sparse, density>,
	    boost::mpl::pair<dense, density>,
	    boost::mpl::pair<banded, layout>,
	    boost::mpl::pair<compressed, layout>,
	    boost::mpl::pair<coordinate, layout>,
	    boost::mpl::pair<ellpack, layout>,
	    boost::mpl::pair<morton, layout>,
	    boost::mpl::pair<symmetric, symmetry>,
	    boost::mpl::pair<anti_symmetric, symmetry>,
	    boost::mpl::pair<self_adjoint, symmetry>,
	    boost::mpl::pair<row_major, orientation>,
	    boost::mpl::pair<col_major, orientation>,
	    boost::mpl::pair<on_stack, location>,
	    boost::mpl::pair<on_heap, location>
	>                            kind_map_common;

	// Find the kind of a parameter
	// templated parameters are handled by specialization
	// Spezialization can be extended in including headers
	template <typename KindMap, typename Kind>
	struct find_kind
	{
	    typedef typename boost::mpl::at<KindMap, Kind>::type type;
	    static_assert( !boost::is_same<type, boost::mpl::void_>::value,
	    		  "The type you providing (the second parameter of find_kind) is not used in the type generator.");
	};
	
	template <typename KindMap, std::size_t ...Values>
	struct find_kind<KindMap, dim<Values...> >
	{
	    typedef dimensionality type;
	};

	template <typename KindMap, std::size_t ...Values>
	struct find_kind<KindMap, mask<Values...> >
	{
	    typedef masking type;
	};

	template <typename KindMap, typename Index>
	struct find_kind<KindMap, as_size_type<Index> >
	{
	    typedef size_type type;
	};
	

	// Error message that this kind of parameter is already set
	template <typename Kind>
	struct error_message_common
	{
	    static_assert( !boost::is_same<Kind, density>::value, "Density (dense, sparse) is set twice.");
	    static_assert( !boost::is_same<Kind, layout>::value, "Layout (banded, compressed, ...) is set twice.");
	    static_assert( !boost::is_same<Kind, symmetry>::value, "Symmetry is set twice.");
	    static_assert( !boost::is_same<Kind, orientation>::value, "Orientation (row_major, col_major) is set twice.");
	    static_assert( !boost::is_same<Kind, location>::value, "Location (on_heap, on_stack) is set twice.");
	    typedef int type;
	};
	
	template <typename TypePara>
	struct replace_defaults_common
	{
	    typedef typename boost::mpl::at<TypePara, dimensionality>::type         dim1;
	    typedef typename boost::mpl::at<TypePara, location>::type               loc1;

	    static_assert( !boost::is_same<loc1, on_stack>::value || boost::is_same<dim1, none>::value,
			  "Types to be stored on stack must provide static size.");

	    // if location not set take stack as default when dim is given otherwise heap
	    typedef typename boost::mpl::if_<
		boost::is_same<loc1, none>,
		typename boost::mpl::if_<boost::is_same<dim1, none>,
					 on_heap,
					 on_stack>::type,
		loc1>::type                                                         loc2;		      

	    typedef typename boost::mpl::at<TypePara, size_type>::type              size1;
	    typedef typename boost::mpl::if_<boost::is_same<size1, none>,
					     as_size_type<std::size_t>,
					     size1>::type                           size2;
	    
	    typedef boost::mpl::map<
		boost::mpl::pair<density, typename boost::mpl::at<TypePara, density>::type>, 
		boost::mpl::pair<layout, typename boost::mpl::at<TypePara, layout>::type>, 
		boost::mpl::pair<symmetry, typename boost::mpl::at<TypePara, symmetry>::type>, 
		boost::mpl::pair<orientation, typename boost::mpl::at<TypePara, orientation>::type>, 
		boost::mpl::pair<dimensionality, dim1>, 
		boost::mpl::pair<size_type, size2>, 
		boost::mpl::pair<location, loc2>,
		boost::mpl::pair<masking, typename boost::mpl::at<TypePara, masking>::type>
	    >                            type;       	
	};


    } // namespace type_para

} // namespace mtl

#endif // defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#endif // MTL_TYPE_PARAMETER_COMMON_INCLUDE
