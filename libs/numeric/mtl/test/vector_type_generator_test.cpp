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

#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/mtl/mtl.hpp>


int main(int, char**)
{
    using namespace mtl;
    using mtl::io::tout;

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_STATICASSERT) && defined(MTL_WITH_TEMPLATE_ALIAS)

    vector<double> v1;

    tout << "Type of v1 is " << typeid(v1).name() << '\n';
    static_assert(!traits::is_row_major<decltype(v1)>::value, "Vector must be column-major!");

    vector<double, column_major> v2;

    tout << "Type of v2 is " << typeid(v2).name() << '\n';
    static_assert(!traits::is_row_major<decltype(v2)>::value, "Vector must be column-major!");

    vector<double, row_major> v3;

    tout << "Type of v3 is " << typeid(v3).name() << '\n';
    static_assert(traits::is_row_major<decltype(v3)>::value, "Vector must be row-major!");

    vector<double, as_size_type<unsigned> > v4;

    tout << "Type of v4 is " << typeid(v4).name() << '\n';
    static_assert(boost::is_same<mtl::Collection<decltype(v4)>::size_type, unsigned>::value, 
		  "Size type must be unsigned!");

    vector<double, dim<3> > v5;

    tout << "Type of v5 is " << typeid(v5).name() << '\n';
    static_assert(traits::is_static<decltype(v5)>::value, "Vector must have static size!");
    

    vector<double, sparse> v6;

    tout << "Type of v6 is " << typeid(v6).name() << '\n';
    static_assert(traits::is_sparse<decltype(v6)>::value, "Vector must be sparse!");


#else
    tout << "Test deactivated due to missing C++11 features.\n";
#endif

    return 0;
}
