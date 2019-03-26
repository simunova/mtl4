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

// With contributions from Cornelius Steinhardt

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/secular.hpp>
#include <boost/numeric/mtl/operation/sort.hpp>

using namespace std;


int main(int , char**)
{
    using namespace mtl;
    typedef dense_vector<double>  Vector;
    Vector                    lambda(2, 0.0), z(2), d(2);
    z[0]=1; z[1]=1;
    d[0]=-5; d[1]=-1;
   
    //lambda= mtl::secular_f<Vector>(lambda, z, d, 5.0).f(3.0);

    mtl::vec::secular_f<Vector>   ss(z, d, 5.0);
    std::cout<<"lambda  =" << ss.f(3.0) <<"\n";
    std::cout<<"lambda  =" << ss.f(0.0) <<"\n";
    std::cout<<"lambda  =" << ss.f(-3.0) <<"\n";
    std::cout<<"lambda  =" << ss.f(13.0) <<"\n";
    std::cout<<"lambda  =" << ss.f(113.0) <<"\n";

    std::cout<<"lambda  =" << ss.grad_f(13.0) <<"\n";
    std::cout<<"lambda  =" << ss.grad_f(113.0) <<"\n";
    std::cout<<"roots  =" << secular(z, d, 5.0) <<"\n";
    //std::cout<<"lambda  =" << lambda <<"\n";

    Vector x(5, 0.0);
    for(int i = 0; i < 5; i++)
	x[i]=5-i;
    x[1]=1;
    std::cout<< "\n x=" << x << "\n";
    sort(x);
    std::cout<< "x=" << x << "\n";
    MTL_THROW_IF(x[0] != 1.0, mtl::runtime_error("Error in sorting."));

    return 0;
}



