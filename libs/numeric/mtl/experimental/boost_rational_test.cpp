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
#include <boost/rational.hpp>

int main(int argc, char* argv[]) 
{
#if 0
    boost::rational<long> test(1,2);
    test += boost::rational<long>(1,3);
    std::cout << test;

    mtl::dense_vector< boost::rational<long> > testVec(2);
    testVec(0) = boost::rational<long>(1,3);
    testVec[1] = boost::rational<long>(1,4);

    std::cout << testVec(0) << "\n";
    std::cout << testVec[1] << "\n";
    std::cout << testVec << "\n";
#endif
}
