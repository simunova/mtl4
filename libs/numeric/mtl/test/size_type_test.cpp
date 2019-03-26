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

// Tests only compilability with pedantic warnings, run-time values aren't tested

#include <cstddef>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

template <typename Size>
void test()
{
    io::tout << "\n\nSize type is " << typeid(Size).name() << '\n';
    typedef mat::parameters<row_major, mtl::index::c_index, non_fixed::dimensions, false, Size> para;
    compressed2D<double, para>   A;
    laplacian_setup(A, 2, 3);
    
    io::tout << "A is\n" << A << "\nsize of index is " << sizeof(A.ref_minor()[0]) << '\n';
}

int main()
{
    test<boost::uint_least16_t>();
    test<boost::uint_least32_t>();
    test<std::size_t>();

    test<boost::int_least16_t>();
    test<boost::int_least32_t>();
    test<std::ptrdiff_t>();

    return 0;
}
