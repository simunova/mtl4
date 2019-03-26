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

// Regression test for problem reported by Garth Wells

#include <boost/numeric/mtl/mtl.hpp>

namespace mtl {
    using mtl::num_rows; // has no effect with the friend definition
}


int main(int, char**)
{
    mtl::dense_vector<double> x(10);
    std::size_t size1 = num_rows(x);
    std::size_t size2 = mtl::num_rows(x);      // does not compile with friend definition
    std::size_t size3 = mtl::num_rows(x);      // does not compile with friend definition either

    std::cout << size1 + size2 + size3 << "\n";

    return 0; 
}
