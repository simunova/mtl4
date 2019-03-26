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
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// To improve our coverage ;-)

template <typename Exception>
void test(Exception)
{
    try {
	throw Exception();
    }
    catch (Exception e) {}
}



int main()
{
    using namespace mtl;
    test(index_out_of_range());
    test(mtl::range_error());
    test(mtl::domain_error());
    test(incompatible_size());
    test(need_nonempty());
    test(change_static_size());
    test(argument_result_conflict());
    test(incompatible_size());
    test(matrix_not_square());
    test(matrix_too_small());
    test(matrix_singular());
    test(missing_diagonal());
    test(access_during_insertion());
    test(unexpected_result());
    test(mtl::runtime_error());
    test(mtl::logic_error());
    test(io_error());
    test(file_not_found());

#if 0
    int a= 3, b= 4;
#ifndef MTL_ASSERT_FOR_THROW
#  define MTL_ASSERT_FOR_THROW
#endif
    MTL_DEBUG_THROW_IF(a != b, unexpected_result()); // to cover assertion in exceptions
#endif

    return 0;
}
