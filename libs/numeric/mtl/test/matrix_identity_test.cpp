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
    
    vt v(60);
    iota(v);
    tout << "v is " << v << endl;

    mtl::mat::identity2D I(60);
    vt w1(I * v);
    tout << "I * v is " << w1 << endl;

    if (one_norm(vt(w1 - v)) > 0.001) throw "Wrong result with square identity";

    w1-= I * v;
    tout << "w1-= I * v is " << w1 << endl;
    if (one_norm(w1) > 0.001) throw "Wrong result";

    w1+= I * v;
    tout << "w1+= I * v is " << w1 << endl;
    if (one_norm(vt(w1 - v)) > 0.001) throw "Wrong result with square identity"; 

    vt w2( w1 - I * v );
    double alpha;
    (lazy(w2)= I * v) || (lazy(alpha)= lazy_dot(w2, v));

    vt w3(30), w4(90);
    mtl::mat::identity2D I3(30, 60), I4(90, 60);
    
    w3= I3 * v;
    tout << "I3 * v is " << w3 << endl;
    if (one_norm(vt(w3 - v[mtl::irange(30)])) > 0.001) throw "Wrong result with broad identity";
    
    w4= I4 * v;
    tout << "I4 * v is " << w4 << endl;
    if (one_norm(vt(w4[mtl::irange(60)] - v)) > 0.001) throw "Wrong result with long identity";
    if (one_norm(w4[mtl::irange(60, 90)]) > 0.001) throw "Wrong result with long identity";

#if 0
    mtl::dense2D<double> A(60,60);
    A=2;
    A+= A*I;
    std::cout<< "A=\n" << A << "\n";
#endif
    return 0;
}
