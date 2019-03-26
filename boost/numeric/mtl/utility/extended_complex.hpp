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

#ifndef MTL_TRAITS_EXTENDED_COMPLEX_INCLUDE
#define MTL_TRAITS_EXTENDED_COMPLEX_INCLUDE

#include <complex>
#include <boost/numeric/mtl/utility/different_non_complex.hpp>
#include <boost/numeric/mtl/concept/std_concept.hpp>

namespace mtl { namespace traits {


template <typename T, typename U, bool enable>
struct extended_complex_aux 
{};

template <typename T, typename U>
struct extended_complex_aux<T, U, true>
{
    typedef std::complex<typename mtl::Addable<T, U>::result_type> type;
};


/// Result type of extended complex binary arithmetic 
/** If operation already exist in standard or non-complex are incompatible \p type does not exist.
    extended_complex<T,U>::type is an implicit enable_if and avoids ambiguities.  **/
template <typename T, typename U>
struct extended_complex
  : extended_complex_aux<T, U, different_non_complex<T, U>::value >
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_EXTENDED_COMPLEX_INCLUDE
