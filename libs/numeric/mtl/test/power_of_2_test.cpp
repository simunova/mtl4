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

#include <boost/numeric/meta_math/power_of_2.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

void test_power_of_2() 
{
    using namespace meta_math;

    MTL_THROW_IF(power_of_2<0>::value != 1, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(power_of_2<1>::value != 2, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(power_of_2<2>::value != 4, mtl::runtime_error("wrong value"));
    MTL_THROW_IF(power_of_2<3>::value != 8, mtl::runtime_error("wrong value"));
}


int main() 
{
    test_power_of_2();
    return 0;
}
