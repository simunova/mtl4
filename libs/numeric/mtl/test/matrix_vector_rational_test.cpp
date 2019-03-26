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
typedef mtl::dense_vector<t_Q> t_dVecQ;
typedef mtl::dense2D<t_Q> t_dMatQ;

int main(int, char* [])
{
    //! Test Matrix-Vector operators
    t_dMatQ mQ1(3,3); mtl::mat::diagonal_setup(mQ1,2);
    std::cout << "mQ1:\n" << mQ1 << "\n";
    std::cout << "size(mQ1):\n" << mtl::mat::size(mQ1) << "\n";

    t_dVecQ vQ4(3,1);
    //vQ4 *= mQ1;      // not defined yet -> ticket #254
    std::cout << "vQ4: " << vQ4 << "\n";

    t_dVecQ vQ5( mQ1 * vQ4 );
    std::cout << "vQ5: " << vQ5 << "\n";

    return 0;
}

