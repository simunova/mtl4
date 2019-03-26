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


inline void bench(unsigned n)
{
    typedef mtl::mat::parameters<> mat_para;
    mtl::dense2D<double, mat_para> A(n, n), B(n, n), C(n, n);
    mtl::dense_vector<double> x(n, 1.0), b(n);

    for ( unsigned i= 0; i < n; i++)
	for ( unsigned  j= 0; j < n; j++){
	    A[i][j]= 1;
	    B[i][j]= 2;
	}
    
    boost::timer t;
    C= A * B;
    b= C * x;
    double t1= 1000. * t.elapsed();

    t.restart();
    b=  B * x;
    x= A * b;
    double t2= 1000. * t.elapsed();
    std::cout << n << " " << t1 << " " << t2 <<"\n";
}



int main(int , char**)
{
    bench(1000);



#if 0
    bench(300);
    bench(400);
    bench(500);
    bench(600);
    bench(700);
    bench(800);
    bench(900);
    bench(1000);
    bench(1200);
    bench(1400);
    bench(1600);
    bench(1800);
    bench(2000);
    bench(2200);
    bench(2400);
    bench(2600);
#endif

    return 0;
}
