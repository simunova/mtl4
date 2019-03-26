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

#ifndef MTL_TYPE_PARAMETER_INCLUDE
#define MTL_TYPE_PARAMETER_INCLUDE

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#include <boost/numeric/mtl/utility/type_parameter_local.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/erase_key.hpp>
#include <boost/mpl/insert.hpp>

namespace mtl {

    namespace type_para {

	typedef init_local        init;
	typedef kind_map_local    kind_map;
	
	// Error message that this kind of parameter is already set
	template <typename Kind>
	struct error_message
	  : error_message_local<Kind> {};

	// Print an error message if the value is not none, i.e. the parameter is already set
	template <typename Element>
	struct check
	  : error_message<typename Element::first> {};

	// Spezialize for virgin parameters, nothing to do
	template <typename Key>
	struct check<boost::mpl::pair<Key, none> > { typedef int type; };


	// Set parameter Value by replacing according entry in TypePara
	template <typename TypePara, typename Value>
	struct set_para
	{
	    typedef typename find_kind<kind_map, Value>::type       kind;

	    // find entry in type paras
	    typedef typename boost::mpl::at<TypePara, kind>::type   value;
	    
	    // check whether is not yet set
	    typedef typename check<boost::mpl::pair<kind, value> >::type    check_dummy;

	    // erase empty entry
	    typedef typename boost::mpl::erase_key<TypePara, kind>::type    short_para;
	   
	    // insert new entry
	    typedef typename boost::mpl::insert<short_para, short_para, boost::mpl::pair<kind, Value> >::type type;
	};

	
	template <typename TypePara, typename ...Values>
	struct set_parameters
	{
	    typedef TypePara       type;
	};

	template <typename TypePara, typename FirstValue, typename ...Values>
	struct set_parameters<TypePara, FirstValue, Values...>
	{
	    typedef typename set_para<TypePara, FirstValue>::type              next_parameters;
	    typedef typename set_parameters<next_parameters, Values...>::type  type;
	};

	template <typename TypePara>
	struct replace_defaults
	  : replace_defaults_local<TypePara> {};
	


    } // namespace type_para

    template <typename ...Values>
    struct set_parameters
    {
	typedef typename type_para::set_parameters<type_para::init, Values...>::type parameters_by_user;
	typedef typename type_para::replace_defaults<parameters_by_user>::type       type;
    };


} // namespace mtl

#endif // defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#endif // MTL_TYPE_PARAMETER_INCLUDE
