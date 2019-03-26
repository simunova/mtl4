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

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/io/test_ostream.hpp>

using namespace std;




int main()
{
    using mtl::io::tout;
    using mtl::vec::parameters;

    bool with_errors= false;

    mtl::dense_vector<int>                                u(3), v, ref(3);
    mtl::dense_vector<int, parameters<mtl::row_major> >   ur(3), vr, refr(3);

    u= 1, 2, 3;
    ur= 1, 2, 3;

    v= u * u;   
    tout << "u * u = " << v << endl;

    vr= ur * ur;   
    tout << "ur * ur = " << vr << endl;

    ref= 1, 4, 9;
    refr= 1, 4, 9;

    if (v != ref) {
	cerr << "Wrong result in element-wise product (column vectors)" << endl;
	with_errors= true;
    }

    if (vr != refr) {
	cerr << "Wrong result in element-wise product (row vectors)" << endl;
	with_errors= true;
    }
    
    return with_errors ? 1 : 0;
}
