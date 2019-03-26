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

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/is_static.hpp>


int main(int, char**)
{
    using namespace mtl;

    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3>, true> mat_para;
    typedef mtl::vec::parameters<tag::col_major, mtl::vec::fixed::dimension<3>, true>                vec_para;

    MTL_THROW_IF(( traits::is_static<mtl::non_fixed::dimensions>::value), mtl::runtime_error("Must not be static!"));
    MTL_THROW_IF((!traits::is_static<mtl::fixed::dimensions<1, 2> >::value), mtl::runtime_error("Must be static!"));

    MTL_THROW_IF(( traits::is_static<mtl::vec::non_fixed::dimension>::value), mtl::runtime_error("Must not be static!"));
    MTL_THROW_IF((!traits::is_static<mtl::vec::fixed::dimension<1> >::value), mtl::runtime_error("Must be static!"));

    MTL_THROW_IF(( traits::is_static<mtl::dense2D<float> >::value), mtl::runtime_error("Must not be static!"));
    MTL_THROW_IF((!traits::is_static<mtl::dense2D<float, mat_para> >::value), mtl::runtime_error("Must be static!"));

    MTL_THROW_IF(( traits::is_static<mtl::morton_dense<float, morton_mask> >::value), mtl::runtime_error("Must not be static!"));
    MTL_THROW_IF((!traits::is_static<mtl::morton_dense<float, morton_mask, mat_para> >::value), mtl::runtime_error("Must be static!"));

    MTL_THROW_IF(( traits::is_static<mtl::compressed2D<float> >::value), mtl::runtime_error("Must not be static!"));
    MTL_THROW_IF((!traits::is_static<mtl::compressed2D<float, mat_para> >::value), mtl::runtime_error("Must be static!"));

    MTL_THROW_IF(( traits::is_static<mtl::dense_vector<float> >::value), mtl::runtime_error("Must not be static!"));
    MTL_THROW_IF((!traits::is_static<mtl::dense_vector<float, vec_para> >::value), mtl::runtime_error("Must be static!"));

    return 0;
}
