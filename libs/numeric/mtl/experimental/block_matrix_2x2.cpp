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

    using mblock= matrix<double, dim<2, 2> >;
    using vblock= mtl::vector<double, dim<2> >;

    using mtype= matrix<mblock, sparse>;
    using vtype= mtl::vector<vblock>;

    mtype A(2, 3);
    {
	mat::inserter<mtype> ins(A);
	mblock mb;
	mb[0][0]= 1; mb[0][1]= 2; 
	mb[1][0]= 3; mb[1][1]= 4;

	ins[0][0] << mb;

	mb[1][1]= 5;
	ins[1][1] << mb;

	// ins(0, 0) << mblock{{1, 2}, {3, 4}};
	// ins[1][1] << mblock{{3, 4}, {5, 6}};
    }
    // cout << "A = " << A;

    vtype x{vblock{1, 3}, vblock{1, 2}, vblock{9, 3}};
    cout << "x = " << x << endl;

    vtype y( A * x );
    cout << "y= " << y << endl;

    return 0;
}
 
