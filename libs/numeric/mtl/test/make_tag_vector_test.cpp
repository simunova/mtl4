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

int main()
{
    using namespace mtl;
    size_t  N= 18;

    dense_vector<bool> v1( make_tag_vector(N, mtl::srange(0, imax, 2)) );
    std::cout << "Strided vector v1 is " << v1 << '\n';

    MTL_THROW_IF(size(v1) != N, runtime_error("Wrong size"));
    for (size_t i= 0; i < N; i++)
	MTL_THROW_IF(v1[i] != (i % 2 == 0), runtime_error("Wrong value for strided tags"));

    dense_vector<bool> v2( make_tag_vector(N, mtl::irange(7, 11)) );
    std::cout << "Strided vector v2 is " << v2 << '\n';

    MTL_THROW_IF(size(v2) != N, runtime_error("Wrong size"));
    for (size_t i= 0; i < N; i++)
	MTL_THROW_IF(v2[i] != (i >= 7 && i < 11), runtime_error("Wrong value for ranged tags"));

    dense_vector<bool> v3( make_tag_vector(N, mtl::irange(7, imax)) );
    std::cout << "Strided vector v3 is " << v3 << '\n';

    MTL_THROW_IF(size(v3) != N, runtime_error("Wrong size"));
    for (size_t i= 0; i < N; i++)
	MTL_THROW_IF(v3[i] != (i >= 7), runtime_error("Wrong value for ranged tags"));

    return 0;
}
