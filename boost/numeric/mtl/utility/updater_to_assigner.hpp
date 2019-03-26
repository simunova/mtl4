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

#ifndef MTL_TRAITS_UPDATER_TO_ASSIGNER_INCLUDE
#define MTL_TRAITS_UPDATER_TO_ASSIGNER_INCLUDE

#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/operation/update.hpp>

namespace mtl { namespace traits {

template <typename Updater> 
struct updater_to_assigner
{};

template <typename Element>
struct updater_to_assigner<mtl::operations::update_store<Element> >
{
    typedef mtl::assign::assign_sum type;
};

template <typename Element>
struct updater_to_assigner<mtl::operations::update_plus<Element> >
{
    typedef mtl::assign::plus_sum type;
};

template <typename Element>
struct updater_to_assigner<mtl::operations::update_minus<Element> >
{
    typedef mtl::assign::minus_sum type;
};

}} // namespace mtl::traits

#endif // MTL_TRAITS_UPDATER_TO_ASSIGNER_INCLUDE
