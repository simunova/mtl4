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

#ifndef MTL_UTILITY_SOMETIMES_DATA_INCLUDE
#define MTL_UTILITY_SOMETIMES_DATA_INCLUDE

namespace mtl { namespace utility {

template <bool C, typename T>
struct sometimes_data
{
    sometimes_data(const T& data) : data(data) {}
    T data;
};

template <typename T>
struct sometimes_data<false, T>
{
    sometimes_data(const T&) {}
};

} // namespace utility

using utility::sometimes_data;

} // namespace mtl

#endif // MTL_UTILITY_SOMETIMES_DATA_INCLUDE
