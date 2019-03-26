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
#include <boost/timer.hpp>
#include <cstdlib>


using namespace std;

void test(int n, int s)
{
    mtl::compressed2D<double> A(n, n), B(n, n);
    {
	double val = 3.14;
	mtl::mat::inserter<mtl::compressed2D<double> > ins(A, int(1.2*s)), ins2(B, int(1.2*s));
	for (int i= 0; i < n*s; i++) {
	    ins[rand()%n][rand()%n] << val;
	    ins2[rand()%n][rand()%n] << val;
	}
    }
    boost::timer start;
    if (n < 10) cout << "A is \n" << A << "B is \n" << B;
    A+= B;
    if (n < 10) cout << "A is \n" << A;
    cout << "Addition took " << start.elapsed() << "s\n";
}

 
int main(int argc, char* argv[])
{
    test(atoi(argv[1]), atoi(argv[2]));

    return 0; 
}
