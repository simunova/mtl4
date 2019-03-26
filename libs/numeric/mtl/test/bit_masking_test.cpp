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

#include <cstdio>
#include <iostream>
#include <boost/numeric/mtl/recursion/bit_masking.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>


using namespace std;  


template <unsigned long Mask>
void test()
{
    printf("Mask %lx, 32x32 base is row-major %i, is column-major %i, shark 2 row-major %i"
	   ", 4x4 row-major %i, column-major %i\n",
	   Mask, mtl::is_32_base_case_row_major<Mask>::value,
	   mtl::is_k_power_base_case_col_major<5, Mask>::value,
	   mtl::is_k_power_base_case_row_major_t_shark<5, 1, Mask>::value,
	   mtl::is_k_power_base_case_row_major<2, Mask>::value,
	   mtl::is_k_power_base_case_col_major<2, Mask>::value);
}

template <unsigned long Mask1, unsigned long Mask2>
void check_same_mask()
{
    printf("Mask1 %lx, Mask2 %lx\n", Mask1, Mask2);
    MTL_THROW_IF(Mask1 != Mask2, mtl::runtime_error("Different masks\n"));
}

#if 0 // Do I use this???
template <bool is_4, unsigned long long s, unsigned long long l>
struct mm
{
    static const unsigned long value= (const unsigned long) l;
};

template <unsigned long long s, unsigned long long l>
struct mm<true, s, l>
{
    static const unsigned long value= (const unsigned long) s;
};


template <unsigned long long s, unsigned long long l>
struct mask
{
    static const unsigned long value= mm<sizeof(unsigned long) == 4, s, l>::value;
};
#endif


int main(int, char**)
{
    using mtl::row_major; using mtl::col_major; using mtl::generate_mask;

    const unsigned long z= 0, morton= (z-1) / 3, morton_z= ~morton, doppled_4_row= morton + 7,
                        doppled_4_col= morton - 2, doppled_32_row= morton + 651, 
	                doppled_32_col= morton - 310, doppled_32_row_shark_2= morton + 620;

    test<morton>();
    test<morton_z>();
    test<doppled_4_row>();
    test<doppled_4_col>();
    test<doppled_32_row>();
    test<doppled_32_col>();
    test<doppled_32_row_shark_2>();

    const unsigned long morton_gen= generate_mask<true, 0, row_major, 0>::value,
	morton_z_gen= generate_mask<false, 0, row_major, 0>::value,
	doppled_4_row_gen= generate_mask<true, 2, row_major, 0>::value,
	doppled_4_col_gen= generate_mask<true, 2, col_major, 0>::value,
	doppled_32_row_gen= generate_mask<true, 5, row_major, 0>::value,
	doppled_32_col_gen= generate_mask<true, 5, col_major, 0>::value,
	doppled_32_row_shark_2_gen= generate_mask<true, 5, row_major, 1>::value;

    check_same_mask<morton, morton_gen>();
    check_same_mask<morton_z, morton_z_gen>();
    check_same_mask<doppled_4_row, doppled_4_row_gen>();
    check_same_mask<doppled_4_col, doppled_4_col_gen>();
    check_same_mask<doppled_32_row, doppled_32_row_gen>();
    check_same_mask<doppled_32_col, doppled_32_col_gen>();
    check_same_mask<doppled_32_row_shark_2, doppled_32_row_shark_2_gen>();

    return 0;
}
