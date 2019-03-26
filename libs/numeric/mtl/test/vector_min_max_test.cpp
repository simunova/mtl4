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
#include <cmath>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/operation/min.hpp>
#include <boost/numeric/mtl/operation/max.hpp>



using namespace std;  
    

template <typename Vector>
void test(Vector& v, const char* name)
{
    typedef typename mtl::Collection<Vector>::value_type value_type;
    using mtl::min; using mtl::max; 

    for (unsigned i= 0; i < size(v); i++)
	v[i]= value_type(double(i+1) * pow(-1.0, int(i))); // Amb. in MSVC 

    std::cout << "\n" << name << "  --- v = " << v; std::cout.flush();

    std::cout << "min(v) = " << min(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(min(v) != -4.0, mtl::runtime_error("min wrong"));

    std::cout << "min<4>(v) = " << min<4>(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(min<4>(v) != -4.0, mtl::runtime_error("min<4> wrong"));

    std::cout << "max(v) = " << max(v) << "\n"; std::cout.flush();
    MTL_THROW_IF(max(v) != 5.0, mtl::runtime_error("max wrong"));
}
 

int main(int, char**)
{
    mtl::dense_vector<float>   u(5);
    mtl::dense_vector<double>  x(5);

    std::cout << "Testing vector operations\n";

    test(u, "test float");
    test(x, "test double");

    mtl::dense_vector<float, mtl::vec::parameters<mtl::row_major> >   ur(5);
    test(ur, "test float in row vector");

    return 0;
}
 














