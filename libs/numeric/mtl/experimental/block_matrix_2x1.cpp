// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschrÃ¤nkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;

int main(int argc, char** argv) 
{
    using namespace mtl;

    using cblock= mtl::vector<double, dim<2> >;
    using rblock= mtl::vector<double, dim<2>, row_major >; 

    using mtype= matrix<rblock, sparse>;
    using vtype= mtl::vector<cblock>;

    mtype A(2, 3);
    {
	mat::inserter<mtype> ins(A);
	
	ins[0][0] << rblock{1, 3};
	ins[1][1] << rblock{4, 9};
    }
    cout << "A = \n" << A;

    vtype x{cblock{1, 3}, cblock{1, 2}, cblock{9, 3}};
    cout << "x = " << x << endl;

    mtl::vector<double> y( A * x );
    cout << "y= " << y << endl;

    return 0;
}
 
