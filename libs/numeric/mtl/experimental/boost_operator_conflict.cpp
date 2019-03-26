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

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/operators.hpp>

class test : boost::addable<test>
{
  public:
    test& operator+=(const test &op) { }
};

int main()
{
    mtl::dense_vector<test> v(7);
    for(int i = 0; i < mtl::size(v); i++)
	;
    return 0;
}
