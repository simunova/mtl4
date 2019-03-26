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

#ifndef MTL_TRAITS_VIEW_CODE_INCLUDE
#define MTL_TRAITS_VIEW_CODE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/is_what.hpp>

namespace mtl { namespace traits {

/// Represent views by code (quite advanced feature for internal use)
/** 1-bit: constant
    2-bit: conjugated
    4-bit: transposed
    -> hermitian_view of constant matrix is 7. **/
template <typename T>
struct view_code
{
    MTL_STATIC_ASSERT(is_matrix<T>::value, "Currently only matrices are supported.");
    static const unsigned value= 0; ///< Default is zero
};

template <typename T>
struct view_code<const T>
{
    static const unsigned value= view_code<T>::value | 1; ///< Toggle constancy bit
};

template <typename Matrix>
struct view_code< mat::conj_view<Matrix> >
{
    static const unsigned value= view_code<Matrix>::value ^ 2; ///< Toggle conjugation bit
};

template <typename Matrix> 
struct view_code<mtl::mat::transposed_view<Matrix> >
{
    static const unsigned value= view_code<Matrix>::value ^ 4; ///< Toggle transposition bit
};

template <typename Matrix> 
struct view_code<mtl::mat::hermitian_view<Matrix> >
{
    static const unsigned value= view_code<Matrix>::value ^ 6; ///< Toggle transposition and conjugation bit
};

// add vector stuff

template <unsigned Value>
struct view_normalize_const
{
	static const unsigned tmp2 = Value == 0 || Value == 4 ? Value | 1 : Value; // if matrix ref or transposed, make it const
	static const unsigned value = (tmp2 & 3) == 3 ? tmp2 ^ 1 : tmp2;           // for conj turn off const
};

template <typename ViewCode>
struct view_add_const
{
    static const unsigned value= ViewCode::value | 1;
};

template <typename ViewCode>
struct view_remove_const
{
    static const unsigned value= ViewCode::value & ~1;
};

template <typename ViewCode>
struct view_toggle_conj
{
    static const unsigned value= view_normalize_const<ViewCode::value ^ 2>::value;
};

template <typename ViewCode>
struct view_toggle_trans
{
    static const unsigned value= ViewCode::value ^ 4;
};

template <typename ViewCode>
struct view_toggle_hermitian
{
    static const unsigned value= view_normalize_const<ViewCode::value ^ 6>::value;
};

}} // namespace mtl::traits

#endif // MTL_TRAITS_VIEW_CODE_INCLUDE
