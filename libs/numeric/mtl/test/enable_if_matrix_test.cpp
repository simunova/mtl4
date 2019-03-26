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
#include <string>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>


using namespace std;

enum tag_type { mtag, vtag, stag };


namespace mtl { namespace mat {

    template <typename A1, typename A2>
    typename mtl::traits::enable_if_matrix<A1>::type
    inline f(const A1&, const A2&, const char* name, tag_type tag)
    {
	cout << name << " in namespace mtl::matrix\n";
	MTL_THROW_IF(tag != mtag, mtl::runtime_error("Function shouldn't be called in namespace matrix!"));
    }
}}

namespace mtl { namespace vec {

    template <typename A1, typename A2>
    typename mtl::traits::enable_if_vector<A1>::type
    inline f(const A1&, const A2&, const char* name, tag_type tag)
    {
	cout << name << " in namespace mtl::vector\n";
	MTL_THROW_IF(tag != vtag, mtl::runtime_error("Function shouldn't be called in namespace vector!"));
    }
}}

namespace mtl { namespace scalar {
    template <typename A1, typename A2>
    typename mtl::traits::enable_if_scalar<A1>::type
    inline f(const A1&, const A2&, const char* name, tag_type tag)
    {
	cout << name << " in namespace mtl::vector\n";
	MTL_THROW_IF(tag != stag, mtl::runtime_error("Function shouldn't be called in namespace scalar!"));
    }
}}

namespace mtl {
    using vec::f;
    using mat::f;
    using scalar::f;
}

string inline is(bool b)
{
    return b ? string("is a ") : string("is not a ");
}

int main(int , char**)
{
    using namespace mtl; using namespace mtl::traits;

    typedef mtl::dense2D<double>        m_t;
    typedef mtl::dense_vector<double>   v_t;

    mtl::dense2D<double>       A;
    mtl::dense_vector<double>  v;
    int                        x;

    cout << "A " << is(is_matrix<m_t>::value) << "matrix, "
	 << is(is_vector<m_t>::value) << "vector, "
	 << is(mtl::traits::is_scalar<m_t>::value) << "scalar.\n";

    cout << "v " << is(is_matrix<v_t>::value) << "matrix, "
	 << is(is_vector<v_t>::value) << "vector, "
	 << is(mtl::traits::is_scalar<v_t>::value) << "scalar.\n";

    cout << "x " << is(is_matrix<int>::value) << "matrix, "
	 << is(is_vector<int>::value) << "vector, "
	 << is(mtl::traits::is_scalar<int>::value) << "scalar.\n";

    f(A, A, "f(A, A)", mtag);
    f(A, v, "f(A, v)", mtag);
    f(A, x, "f(A, x)", mtag);

    f(v, A, "f(v, A)", vtag);
    f(v, v, "f(v, v)", vtag);
    f(v, x, "f(v, x)", vtag);

    f(x, A, "f(x, A)", stag);
    f(x, v, "f(x, v)", stag);
    f(x, x, "f(x, x)", stag);

    return 0;
}
