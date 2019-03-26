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


int main(int , char**)
{
    using namespace mtl;
 
    typedef mtl::vec::parameters<tag::col_major, mtl::vec::fixed::dimension<2> > fvec_para;
    typedef mat::parameters<tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<2, 2> > fmat_para;

    dense2D<double, fmat_para>        A; // dimension not needed here
    dense_vector<double, fvec_para>   v, w, w2; // here neither

    A= 1., 2.,
       3., 4.;
    v= 3., 4.;
    w2= 11., 25.;

    w= A * v;
    
    std::cout << "A * v is " << w << "\n\n";
    MTL_THROW_IF(two_norm(w - w2) > 0.001, runtime_error("Matrix vector product wrong"));

    return 0;
}
