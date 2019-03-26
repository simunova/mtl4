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

// Currently needs -DMTL_DEEP_COPY_CONSTRUCTOR !!!

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace std;
    
    mtl::compressed2D<double> m(3, 3);
    {
	mtl::mat::inserter<mtl::compressed2D<double> > ins(m);
	ins(0, 1) << 2.0; ins(1, 0) << 1.0;
	ins(1, 1) << 4.0; ins(2, 2) << 5.0;
    }

    mtl::dense_vector<double> x(3), y(3);
    for (unsigned i= 0; i < size(x); i++) x[i]= double(i+1);

    y = trans(m) * x;
    cout << y << '\n';

    MTL_THROW_IF(y[0] != 2.0, mtl::runtime_error("y[0] should be 2.0!\n"));

    return 0;
}
