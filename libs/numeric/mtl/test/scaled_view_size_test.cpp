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


using namespace std;




int main(int, char**)
{
    typedef mtl::dense_vector<float>         v_type;
    typedef mtl::vec::scaled_view<float, v_type> s_type;

    v_type v(3, 4.0);
    s_type s(2.0f * v);

    cout << size(s) << "\n";
    cout << num_rows(s) << "\n";
    cout << num_cols(s) << "\n";

    return 0;
}



