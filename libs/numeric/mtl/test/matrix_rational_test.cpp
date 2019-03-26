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
# include "boost/rational.hpp"
# include "boost/range.hpp"

typedef boost::rational<long>  t_Q;
typedef mtl::dense2D<t_Q> t_dMatQ;

int main(int, char**)
{
    //! Test Matrix-Matrix operators
    t_dMatQ mQ1(3,3); mtl::mat::diagonal_setup(mQ1,2);
    t_dMatQ mQ2(3,3);
    t_dMatQ mQ3(3,2);
    t_dMatQ mQ4(2,2);
    t_dMatQ mQ5;

    std::cout << "mQ1:\n" << mQ1 << "\n";
    std::cout << "size(mQ1):\n" << mtl::mat::size(mQ1) << "\n";
 
    mQ2 *= mQ1;
    mQ2 += mQ1;
    mQ4 = trans(mQ3)*mQ3;

    mQ5 = t_dMatQ( trans(mQ3)*mQ3 );
    mQ5 = t_dMatQ( mQ4+mQ4 );

    return 0;
}
