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
#include <boost/numeric/mtl/io/write_ast.hpp>


using namespace std;

int main() 
{
    mtl::dense_vector<float> u(3, 1), v(3, 2), w(3, 3), x(3, 4), y;

    y= 4.5 * u + v * trans(w) * x;
    std::cout << "y = " << y << '\n';
 
    mtl::io::write_ast(y+= 4.5 * u + v * trans(w) * x, "ast_example.dot");

    return 0;
}
