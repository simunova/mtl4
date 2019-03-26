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
//
// Written by Jiahu Deng and Peter Gottschling

#include <iostream>
#include <cstdio>
#include <boost/numeric/mtl/detail/dilated_int.hpp>


using namespace std;

struct morton_exception {};

template <typename T>
void test_inc(string s, T dil, typename T::value_type exp_increment) 
{ 
    printf("%s %x, bit mask is %x, negated mask is %x\n", s.c_str(), dil.i, dil.bit_mask, dil.anti_mask);
    if ((++dil).i != exp_increment) throw morton_exception();
 
    for (typename T::value_type i = 1; i < 6; ++i) {
	T dilated(i);
	if (dilated.undilate() != i) throw morton_exception();

	cout << "dilated(" << i << ") = " << dilated <<  "\n";
	if (dil.i != dilated.i) throw morton_exception();
	// check both pre and post increment
	i & 1 ? ++dil : dil++;
    } 
    
    printf("dilated_zero = %x, dilated_one = %x\n", dil.dilated_zero, dil.dilated_one);

    T dec_dil(6);
    for (typename T::value_type i = 6; i > 0; --i) {
	cout << "dilated(" << i << ") = " << dec_dil <<  "\n";
	if (dec_dil != T(i)) throw morton_exception();

	// check both pre and post decrement
	i & 1 ? --dec_dil : dec_dil--;
    }	
}

template <typename T>
void test_plus1(T)
{
    T a(1), b(3), c(4);
    cout << "a = " << a << ", b = " << b << ", a + b = " << a + b << ", c = " << c << "\n";
    if (a + b != c) throw morton_exception();

    cout << "c - b = " << c - b << "\n";
    if (c - b != a) throw morton_exception();
}  
 
template <typename T> 
void test_plus2(T)
{
    T a(22), b(33), c(55);
    cout << "a = " << a << ", b = " << b << ", a + b = " << a + b << ", c = " << c << "\n";
    if (a + b != c) throw morton_exception();

    cout << "c - b = " << c - b << "\n";
    if (c - b != a) throw morton_exception();
}    

template <typename T>
void test_dilated(string s, T dil, typename T::value_type exp_increment) 
{
    test_inc(s, dil, exp_increment);
    test_plus1(dil); test_plus2(dil);
}


int main(int , char**)
{    
    using namespace mtl;

    dilated_int<unsigned, dilated::odd_bits<unsigned>::value, true>     dil1;
    test_dilated( "Odd normalized", dil1, 2);

    dilated_int<unsigned, dilated::even_bits<unsigned>::value, true>     dil2;
    test_dilated( "Even normalized", dil2, 1);

    dilated_int<unsigned, dilated::odd_bits<unsigned>::value, false>     dil3;
    test_dilated( "Odd anti-normalized", dil3, 0x55555557);

    dilated_int<unsigned, dilated::even_bits<unsigned>::value, false>     dil4;
    test_dilated( "Even anti-normalized", dil4, 0xaaaaaaab);

    return 0;
}

