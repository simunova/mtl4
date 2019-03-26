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

#ifndef MTL_RECURSION_UTILITIES_INCLUDE
#define MTL_RECURSION_UTILITIES_INCLUDE

#include <limits>
#include <cmath>

namespace mtl { namespace recursion {


// Splits a number into a next-smallest power of 2 and rest
std::size_t inline first_part(std::size_t n)
{
    if (n == 0) return 0;

    std::size_t i= std::numeric_limits<std::size_t>::max()/2 + 1;

    while(i >= n) i>>= 1;
    return i;
}


// The remainder of first part
std::size_t inline second_part(std::size_t n)
{
    return n - first_part(n);
}


template <typename Matrix>
std::size_t inline outer_bound(Matrix const& matrix)
{
  std::size_t max_dim=std::max(num_rows(matrix), num_cols(matrix)), bound= 1;
  for (; bound < max_dim;) bound<<= 1;
  return bound;
}


template <typename Integral>
Integral inline least_significant_one_bit(Integral x)
{
	// hardcore bit hacking trick
#ifdef _MSC_VER
# pragma warning( disable : 4146 )
#endif
    return x & -x; 
#ifdef _MSC_VER
# pragma warning( default : 4146 )
#endif
}


template <typename Integral>
bool inline is_power_of_2(Integral x)
{
  return x == least_significant_one_bit(x);
}


}} // namespace mtl::recursion

#endif // MTL_RECURSION_UTILITIES_INCLUDE
