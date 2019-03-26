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

#ifndef ITL_PC_IS_IDENTITY_INCLUDE
#define ITL_PC_IS_IDENTITY_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/itl/itl_fwd.hpp>

namespace itl { namespace pc {

template <typename PC>
bool is_identity(const PC&) 
{ return false; }

template <typename Matrix, typename Value>
bool is_identity(const itl::pc::identity<Matrix, Value>&)
{ return true; }

template <typename PC>
struct static_is_identity
  : boost::mpl::false_
{};

template <typename Matrix, typename Value>
struct static_is_identity<itl::pc::identity<Matrix, Value> >
  : boost::mpl::true_
{};

}} // namespace itl::pc

#endif // ITL_PC_IS_IDENTITY_INCLUDE
