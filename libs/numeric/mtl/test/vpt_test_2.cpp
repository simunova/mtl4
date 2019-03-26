// Software License for PMTL
// 
// Copyright (c) 2010 SimuNova UG, www.simunova.com.
// All rights reserved.
// Author: Peter Gottschling
// 
// This file is part of the Parallel Matrix Template Library
// 
// The details are regulated by the EULA at http://www.simunova.com/en/eula
//                             respectively http://www.simunova.com/de/agb.

#define MTL_VPT_LEVEL 2

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  

void helper_function()
{
    mtl::vampir_trace<9990> tracer; 

    std::cout << "In helper_function <id=0998>, it is " << (tracer.is_traced() ? "" : "not ") << "traced\n";
}

