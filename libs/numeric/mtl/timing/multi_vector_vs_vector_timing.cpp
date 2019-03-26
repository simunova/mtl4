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
#include <boost/timer.hpp>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;




template <typename Vector>
void test(const char* name)
{
    using boost::timer;

    for (size_t n= 4096; n < 6000000; n*= 2) {

	Vector                    x(3*n, 10.0), y(3*n, 4.0), x2(3*n, 10.0), y2(3*n, 0.0);
	mtl::multi_vector<Vector> A(n, 3), B(n, 3), C(n, 3), D(n, 3);
	A= 0.0; B= 10.0; C= 4.0;

	boost::timer timer;
	int reps= 0;
	for (; timer.elapsed() < 3.0; reps++) 
	    for (int i= 0; i < 10000; i++)
		A+= 2.0 * B + 3 * C;
	double t1= timer.elapsed();
	std::cout << n << " " << timer.elapsed() / reps; // << "\n";

	boost::timer timer2;
	for (int r= 0; r < reps ; r++) 
	    for (int i= 0; i < 10000; i++)
		y2+= 2.0 * x + 3 * y;
	double t2= timer2.elapsed();
	std::cout << " " << timer2.elapsed() / reps << " == " << (t1 - t2) / t2 * 100 << "%\n";
    }
}

int main(int, char**)
{
    test<mtl::dense_vector<double> >("dense_vector<double>");

    return 0;
}
