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

using namespace std;

int main()
{
    mtl::dense_vector<double> vec(5,0.0);
    mtl::dense_vector<double>& refVec = vec;

    typedef mtl::tag::iter::all iall;
    typedef mtl::traits::range_generator<iall, mtl::dense_vector<double> >::type Iter;

    for (Iter iter(mtl::begin<iall>(refVec)), iend(mtl::end<iall>(refVec)); iter != iend; ++iter)
	cout << *iter << "\n";

    return 0;
}
