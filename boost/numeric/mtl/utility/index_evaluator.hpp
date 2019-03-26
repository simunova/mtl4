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

#ifndef MTL_TRAITS_INDEX_EVALUATOR_INCLUDE
#define MTL_TRAITS_INDEX_EVALUATOR_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>

namespace mtl { namespace traits {

/// Result type of index evaluation
template <typename T> 
struct index_evaluator {};

template <typename T, typename U, typename Assign>
struct index_evaluator<lazy_assign<T, U, Assign> >
  : boost::lazy_enable_if<is_vector<T>,
			  boost::mpl::if_<is_vector<U>,
					  mtl::vec::vec_vec_aop_expr<T, U, Assign>, 
					  mtl::vec::vec_scal_aop_expr<T, U, Assign>
					  >
			  >
{};

// ... more

template <typename T, typename U>
struct index_evaluator<fused_expr<T, U> >
{
    typedef mtl::vec::fused_index_evaluator<T, U> type;
};



}} // namespace mtl::traits

#endif // MTL_TRAITS_INDEX_EVALUATOR_INCLUDE
