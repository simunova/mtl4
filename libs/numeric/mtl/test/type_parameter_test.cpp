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
#include <typeinfo>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/type_parameter.hpp>


int main(int, char**)
{
    using namespace mtl;
    using namespace std;

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT)
    // pipe the output into 'c++filt -t' or something like that
    typedef set_parameters<>::type t1;
    cout << typeid(t1).name() << '\n';
    cout << typeid(set_parameters<sparse>::type).name() << '\n';
    cout << typeid(set_parameters<sparse, compressed>::type).name() << '\n';
    cout << typeid(set_parameters<sparse, banded>::type).name() << '\n';
    cout << typeid(set_parameters<morton, mask<0x5555> >::type).name() << '\n';

    cout << typeid(set_parameters<compressed, as_size_type<unsigned> >::type).name() << '\n';
    cout << typeid(set_parameters<compressed, column_major>::type).name() << '\n';
    cout << typeid(set_parameters<compressed, column_major, self_adjoint>::type).name() << '\n';
    cout << typeid(set_parameters<compressed, column_major, self_adjoint, on_stack>::type).name() << '\n';
    cout << typeid(set_parameters<column_major, dim<3, 3> >::type).name() << '\n';
		   
    // Uncomment the following to see if the right error messages appear, cannot be automated
    // cout << typeid(set_parameters<sparse, compressed, int>::type).name() << '\n'; // Kind not found
    // cout << typeid(set_parameters<banded, compressed>::type).name() << '\n';      // Two of a kind

#else
    cout << "Test deactivated due to missing C++11 features.\n";	   
#endif

    return 0;
}
