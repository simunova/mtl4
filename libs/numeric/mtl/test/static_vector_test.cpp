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
#include <boost/numeric/mtl/utility/static_vector.hpp>

#if defined(MTL_WITH_STATICASSERT) && defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)
namespace mtl {

    template <std::size_t ...Values>
    using dimv= static_vector<std::size_t, Values...>;
}
#endif

int main()
{
    using namespace std;
    using namespace mtl;

#if defined(MTL_WITH_STATICASSERT) && defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)
    typedef static_vector<int, 3, 4, 9, -2>          v1;
    typedef static_vector<std::size_t, 9, 2, 8, 8>   v2;

    std::cout << v1::get<0>::value << '\n';
    std::cout << v1::get<1>::value << '\n';
    std::cout << v1::get<3>::value << '\n';

    // std::cout << v1::get<4>::value << '\n'; // must raise static_assert

    static_assert(v1::get<0>::value == 3, "Wrong value!");
    static_assert(v1::get<1>::value == 4, "Wrong value!");
    static_assert(v1::get<3>::value == -2, "Wrong value!");

    std::cout << v2::get<1>::value << '\n';
    static_assert(v2::get<1>::value == 2, "Wrong value!");
    
    static_assert(boost::is_same<v2, v2::get<0> >::value, "get<0> should be be idempotent!");

# ifndef _MSC_VER // creates problems on VS
    typedef dimv<7, 9> d;

    std::cout << d::get<1>::value << '\n';
    static_assert(d::get<1>::value == 9, "Wrong value!");
# endif

#else
    std::cout << "Test disabled because your compiler does not support all required C++11 features.\n";
#endif

    return 0;
}
