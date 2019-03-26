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

#ifndef MTL_BASE_CASE_TEST_INCLUDE
#define MTL_BASE_CASE_TEST_INCLUDE

#include <algorithm>

namespace mtl { namespace recursion {

// Minimum of dimensions is less or equal to the reference value
struct min_dim_test
{
    min_dim_test(std::size_t comp) : comp(comp) {}

    template <typename Recursator>
    bool operator() (Recursator const& recursator) const
    {
	return std::min(recursator.get_value().num_rows(), 
			recursator.get_value().num_cols()) 
	       <= comp;
    }

private:
    std::size_t  comp;
};


// Minimum of dimensions is less or equal to the reference value
//   and it can't be split into 2 sub-matrices less or equal the ref value
struct undivisible_min_dim_test
{
    undivisible_min_dim_test(std::size_t comp) : comp(comp) {}

    template <typename Recursator>
    bool operator() (Recursator const& recursator) const
    {
	std::size_t min_dim= std::min(recursator.get_value().num_rows(), 
				      recursator.get_value().num_cols()),
	            max_dim= std::max(recursator.get_value().num_rows(),
				      recursator.get_value().num_cols());

	return min_dim <= comp && 2 * min_dim > max_dim;
    }

private:
    std::size_t  comp;
};


// Maximum of dimensions is less or equal to the reference value
struct max_dim_test
{
    max_dim_test(std::size_t comp) : comp(comp) {}

    template <typename Recursator>
    bool operator() (Recursator const& recursator) const
    {
		return std::max(num_rows(*recursator), num_cols(*recursator)) <= comp;
    }

private:
    std::size_t  comp;
};


// Same with compile-time reference value
template <unsigned long BaseCaseSize>
struct max_dim_test_static
{
    static const unsigned long base_case_size= BaseCaseSize;

    template <typename Recursator>
    bool operator() (Recursator const& recursator) const
    {
	return std::max(recursator.get_value().num_rows(), 
			recursator.get_value().num_cols()) 
	       <= BaseCaseSize;
    }
};


// Upper bound of dimensions in recursator is less or equal to the reference value
struct bound_test
{
    bound_test(std::size_t comp) : comp(comp) {}

    template <typename Recursator>
    bool operator() (Recursator const& recursator) const
    {
	return recursator.bound() <= comp;
    }

private:
    std::size_t  comp;
};


// Same with compile-time reference value
template <unsigned long BaseCaseSize>
struct bound_test_static
{
    static const unsigned long base_case_size= BaseCaseSize;

    template <typename Recursator>
    bool operator() (Recursator const& recursator) const
    {
	return recursator.bound() <= base_case_size;
    }
};



}} // namespace mtl::recursion

#endif // MTL_BASE_CASE_TEST_INCLUDE
