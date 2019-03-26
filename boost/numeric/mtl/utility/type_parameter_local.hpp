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

#ifndef MTL_TYPE_PARAMETER_LOCAL_INCLUDE
#define MTL_TYPE_PARAMETER_LOCAL_INCLUDE

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#include <boost/numeric/mtl/utility/type_parameter_common.hpp>

namespace mtl {

    namespace type_para {

	typedef init_common     init_local;
	typedef kind_map_common kind_map_local;

	// Error message that kind is already set
	template <typename Kind>
	struct error_message_local
	  : error_message_common<Kind> {};

	template <typename TypePara>
	struct replace_defaults_local
	  : replace_defaults_common<TypePara> {};


    }

} // namespace mtl

#endif // defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)

#endif // MTL_TYPE_PARAMETER_LOCAL_INCLUDE
