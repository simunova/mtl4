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
#include <cmath>
#include <string>

#include <boost/numeric/mtl/mtl.hpp>
 

using namespace std;  




template <typename Coll>
void test(Coll& x, const char* name)
{
    std::cout << "\n" << name << " is set to:\n";
    mtl::seed<typename mtl::Collection<Coll>::value_type> s;
    random(x, s);
    std::cout << x;
}


int main(int, char**)
{
    using namespace mtl;

    const unsigned size= 3; 

    dense_vector<double>                             dv(size);

    compressed2D<double>                             cdc(size, size);
    compressed2D<std::complex<double> >              ccc(size, size);
    dense2D<double>                                  dc(size, size);
    dense2D<double, mat::parameters<col_major> >  dcc(size, size);
    dense2D<float>                                   fc(size, size);
    morton_dense<double,  morton_mask>               mdc(size, size);
    morton_dense<double, doppled_32_col_mask>        mcc(size, size);

    test(dv, "dense_vector");
    test(cdc, "compressed2D");
    test(ccc, "compressed2D complex");
    test(dc, "dense2D");
    test(dcc, "dense2D col-major");
    test(mdc, "pure Morton");
    test(mcc, "Hybrid col-major");

    return 0;
}
