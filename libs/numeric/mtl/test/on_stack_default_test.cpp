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

int main(int, char**)
{
    using namespace mtl;

    typedef mtl::fixed::dimensions<3, 3> fmdim;
    typedef mtl::non_fixed::dimensions   mdim;

    typedef mtl::vec::fixed::dimension<3>    fvdim;
    typedef mtl::vec::non_fixed::dimension   vdim;

    typedef mat::parameters<tag::row_major, mtl::index::c_index, mdim>   mat_para;
    typedef mat::parameters<tag::row_major, mtl::index::c_index, fmdim>  fmat_para;

    typedef mtl::vec::parameters<tag::col_major, vdim>                   vec_para;
    typedef mtl::vec::parameters<tag::col_major, fvdim>                  fvec_para;

    MTL_THROW_IF( mat_para::on_stack, mtl::runtime_error("Must not be on stack!"));
    MTL_THROW_IF(!fmat_para::on_stack, mtl::runtime_error("Must be on stack!"));

    MTL_THROW_IF( vec_para::on_stack, mtl::runtime_error("Must not be on stack!"));
    MTL_THROW_IF(!fvec_para::on_stack, mtl::runtime_error("Must be on stack!"));

    return 0;
}
