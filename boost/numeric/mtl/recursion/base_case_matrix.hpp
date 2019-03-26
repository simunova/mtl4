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

#ifndef MTL_BASE_CASE_MATRIX_INCLUDE
#define MTL_BASE_CASE_MATRIX_INCLUDE

#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/meta_math/is_power_of_2.hpp>
#include <boost/numeric/meta_math/log_2.hpp>
#include <boost/numeric/mtl/recursion/base_case_test.hpp>
#include <boost/numeric/mtl/recursion/bit_masking.hpp>

namespace mtl { namespace recursion {

template <typename Matrix, typename BaseCaseTest>
struct base_case_matrix
{
    typedef Matrix type;
};

template <typename Elt, unsigned long Mask, typename Parameters, typename BaseCaseTest>
struct base_case_matrix<mtl::mat::morton_dense<Elt, Mask, Parameters>, BaseCaseTest>
{
    MTL_STATIC_ASSERT(meta_math::is_power_of_2<BaseCaseTest::base_case_size>::value, "Static base case size must be power of two");
    static const unsigned long base_case_bits= meta_math::log_2<BaseCaseTest::base_case_size>::value;

    typedef typename boost::mpl::if_<
	is_k_power_base_case_row_major<base_case_bits, Mask>
      , mtl::mat::dense2D<Elt, mat::parameters<row_major> >
      , typename boost::mpl::if_<
	    is_k_power_base_case_col_major<base_case_bits, Mask>
	  , mtl::mat::dense2D<Elt, mat::parameters<col_major> >
          , mtl::mat::morton_dense<Elt, Mask, Parameters>
        >::type
    >::type type;
};


}} // namespace mtl::recursion

#endif // MTL_BASE_CASE_MATRIX_INCLUDE
