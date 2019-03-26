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

#ifndef MTL_FOR_EACH_INCLUDE
#define MTL_FOR_EACH_INCLUDE

namespace mtl { namespace recursion {

// Go recursively down to base case and apply function on it
template <typename Matrix, typename Function, typename BaseCaseTest>
void for_each(mat::recursator<Matrix> const& recursator, Function const& f, BaseCaseTest const& is_base)
{
    if (recursator.is_empty())
	return;

    if (is_base(recursator)) {
	f(*recursator);
	return;
    }

    for_each(recursator.north_west(), f, is_base);
    for_each(recursator.south_west(), f, is_base);
    for_each(recursator.north_east(), f, is_base);
    for_each(recursator.south_east(), f, is_base);
}


// Non-const version
template <typename Matrix, typename Function, typename BaseCaseTest>
void for_each(mat::recursator<Matrix>& recursator, Function const& f, BaseCaseTest const& is_base)
{
    typedef mat::recursator<Matrix> recursator_type;

    if (recursator.is_empty())
	return;

    if (is_base(recursator)) {
	f(recursator.get_value());
	return;
    }

    recursator_type  tmp_nw(recursator.north_west()), tmp_sw(recursator.south_west()),
	             tmp_ne(recursator.north_east()), tmp_se(recursator.south_east());
    for_each(tmp_nw, f, is_base);
    for_each(tmp_sw, f, is_base);
    for_each(tmp_ne, f, is_base);
    for_each(tmp_se, f, is_base);
}


}} // namespace mtl::recursion


#endif // MTL_FOR_EACH_INCLUDE
