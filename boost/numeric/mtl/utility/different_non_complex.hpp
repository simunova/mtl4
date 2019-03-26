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

#ifndef MTL_TRAITS_DIFFERENT_NON_COMPLEX_INCLUDE
#define MTL_TRAITS_DIFFERENT_NON_COMPLEX_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/type_traits/is_scalar.hpp>


namespace mtl { namespace traits {

/// Type trait for different non-complex scalars, i.e. pairs of scalars whose complex equivalents are not supported in binary operations
template <typename T, typename U>
struct different_non_complex
  : boost::mpl::bool_< !boost::is_same<T, U>::value && !boost::is_complex<T>::value && !boost::is_complex<U>::value 
                       && boost::is_scalar<T>::value &&  boost::is_scalar<U>::value>
{};


}} // namespace mtl::traits

#endif // MTL_TRAITS_DIFFERENT_NON_COMPLEX_INCLUDE
