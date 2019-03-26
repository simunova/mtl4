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

#include <typeinfo>
#include <iostream>
#include <cassert>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace std;
    using mtl::lazy; using mtl::io::tout;
    typedef mtl::dense_vector<double> vt;
    
    mtl::compressed2D<double> A0;
    laplacian_setup(A0, 4, 15);
    tout << "A0 is\n" << A0 << endl;

    vt v(60);
    iota(v);
    tout << "v is " << v << endl;

    vt w1(A0 * v);
    tout << "A0 * v is " << w1 << endl;

    mtl::mat::poisson2D_dirichlet A(4, 15);
    vt  w2(60);
    w2= A * v;
    tout << "A * v is " << w2 << endl;

    if (one_norm(vt(w1 - w2)) > 0.001) throw "Wrong result";

    w2+= A * v;
    tout << "w2+= A * v is " << w2 << endl;
    if (one_norm(vt(w1 + w1 - w2)) > 0.001) throw "Wrong result";

    w2-= A * v;
    tout << "w2-= A * v is " << w2 << endl;
    if (one_norm(vt(w1 - w2)) > 0.001) throw "Wrong result";

    vt w3( w2 - A * v );
    double alpha;
    (lazy(w3)= A * v) || (lazy(alpha)= lazy_dot(w3, v));

    return 0;
}
