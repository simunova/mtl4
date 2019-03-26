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

#ifndef MTL_BIT_MASKING_INCLUDE
#define MTL_BIT_MASKING_INCLUDE

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include <boost/numeric/mtl/utility/tag.hpp>

namespace mtl {

/*
     The bit masks are row masks, that mean 1s represent rows and 0s columns.

     Bit masks:

     i-order (cyrillic i):
     ---------------------
  
     binary:     01010101 ... 01
     0x55555555


     z-order:
     --------
     
     binary:     10101010 ... 10
     0xaaaaaaaa

     row major:
     ----------

     with 2^k columns
     binary:    111111111...1000...0
                             ------- k 0s at the end (LSB), all bits before 1s (MSB)

     column major:
     -------------

     with 2^k rows
     binary:    000000000...0111...1
                             ------- k 1s at the end (LSB), all bits before 0s (MSB)

     hybrid (Doppled):
     -----------------

     i-order
     with 2^k by 2^k base case
     row major
     binary     0101....011...10...0
                               ----- k 0s at the end (LSB); means columns
                          ----- k 1s before; means rows
                ---------- i order
     e.g. 32 by 32 base case 0101...01 11111 00000 = 0x555557e0

     column major
     binary     0101....010...01...1
                               ----- k 1s at the end (LSB); means rows
                          ----- k 0s before; means columns
                ---------- i order
     e.g. 32 by 32 base case 0101...01 00000 11111 = 0x5555541f

     Shark-tooth base case:
     ----------------------

     2^t tooth length
     in  2^k by 2^k base case (of course t <= k)
     row-major
     binary     1..1 00..0 1..1
                           ---- t 1s at the end (LSB); means 2^t tooth allong rows
                     ---- k 0s before; means columns
                ---- k-t 1s before; means rows

     column-major
     binary     0..0 11..1 0..0
                           ---- t 0s at the end (LSB); means 2^t tooth allong columns
                     ---- k 1s before; means rows
                ---- k-t 0s before; means columns


*/


// Mask for the last N bits
template <unsigned long N>
struct lsb_mask
{
    static const unsigned long value= (lsb_mask<N-1>::value << 1) | 1;
};


template <>
struct lsb_mask<0>
{
    static const unsigned long value= 0;
};


/// Last N bits of Value
template <unsigned long N, unsigned long Value>
struct lsb_bits
{
    static const unsigned long value= lsb_mask<N>::value & Value;
};


/// Compares two masks
template <unsigned long Mask1, unsigned long Mask2>
struct same_mask
{
    static const bool value= false;
};

template <unsigned long Mask>
struct same_mask<Mask, Mask>
{
    static const bool value= true;
};


/// Row-major mask for 2^K by 2^K base case
template <unsigned long K>
struct row_major_mask
{
    static const unsigned long value= lsb_mask<K>::value << K;
};


/// Column-major mask for 2^K by 2^K base case
template <unsigned long K>
struct col_major_mask
    : public lsb_mask<K>
{};


/// Checks whether 2^K by 2^K base case of hybric matrix, defined by Mask, is a row-major matrix
template <unsigned long K, unsigned long Mask>
struct is_k_power_base_case_row_major
{
    static const bool value= same_mask<lsb_bits<2*K, Mask>::value, row_major_mask<K>::value>::value;
    // typedef 
};


/// Checks whether 2^K by 2^K base case of hybric matrix, defined by Mask, is a column-major matrix
template <unsigned long K, unsigned long Mask>
struct is_k_power_base_case_col_major
{
    static const bool value= same_mask<lsb_bits<2*K, Mask>::value, col_major_mask<K>::value>::value;
};


/// Checks whether 32x32 base case of hybric matrix, defined by Mask, is a row-major matrix
template <unsigned long Mask>
struct is_32_base_case_row_major
    : public is_k_power_base_case_row_major<5, Mask>
{};


/// Checks whether 32x32 base case of hybric matrix, defined by Mask, is a col-major matrix
template <unsigned long Mask>
struct is_32_base_case_col_major
    : public is_k_power_base_case_col_major<5, Mask>
{};


/// Row-major mask for 2^K by 2^K base case with 2^T shark teeth
template <unsigned long K, unsigned long T>
struct row_major_shark_mask
{
    static const unsigned long value= (lsb_mask<K-T>::value << (K+T)) | lsb_mask<T>::value;
};


/// Row-major mask for 2^K by 2^K base case with 2^T shark teeth
template <unsigned long K, unsigned long T>
struct col_major_shark_mask
{
    static const unsigned long value= lsb_mask<K>::value << T;
};


/** Checks whether 2^K by 2^K base case of hybric matrix, defined by Mask,
    is a row-major matrix shark-tooth with 2^T tooth length
**/
template <unsigned long K, unsigned long T, unsigned long Mask>
struct is_k_power_base_case_row_major_t_shark
{
    static const bool value= same_mask<lsb_bits<2*K, Mask>::value, row_major_shark_mask<K, T>::value>::value;
};


/** Checks whether 2^K by 2^K base case of hybric matrix, defined by Mask,
    is a col-major matrix shark-tooth with 2^T tooth length
**/
template <unsigned long K, unsigned long T, unsigned long Mask>
struct is_k_power_base_case_col_major_t_shark
{
    static const bool value= same_mask<lsb_bits<2*K, Mask>::value, col_major_shark_mask<K, T>::value>::value;
};

  // e-order
/// N-order mask of N bits
template <unsigned long N>
struct i_order_mask
{
    // Check if N is even !!!
    static const unsigned long value= (i_order_mask<N-2>::value << 2) | 1;
};

template<> struct i_order_mask<0> : public lsb_mask<0> {};  // set to 0


/// Z-order mask of N bits
template <unsigned long N>
struct z_order_mask
{
    // Check if N is even !!!
    static const unsigned long value= (z_order_mask<N-2>::value << 2) | 2;
};

template<> struct z_order_mask<0> : public lsb_mask<0> {};  // set to 0


/** Generate arbitrary hybrid mask.
    \param IOrder if true then i-order otherwise z-order
    \param K      2^K by 2^K base case 
    \param Orientation  mtl::row_major or mtl::col_major
    \param T      2^T tooth length
**/
template <bool IOrder, unsigned long K, typename Orientation, unsigned long T>
class generate_mask
{
    static const unsigned long rec_size= 8 * sizeof(unsigned long) - 2 * K,
	rec_part= (IOrder ? i_order_mask<rec_size>::value : z_order_mask<rec_size>::value) << 2*K;
    typedef typename boost::mpl::if_<
	boost::is_same<Orientation, row_major>
      , row_major_shark_mask<K, T>
      , col_major_shark_mask<K, T>
    >::type base_part_type;
public:
    static const unsigned long value= rec_part | base_part_type::value;
};


} // namespace mtl

#endif // MTL_BIT_MASKING_INCLUDE
