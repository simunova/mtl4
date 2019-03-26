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

#ifndef META_MATH_LOOP2_INCLUDE
#define META_MATH_LOOP2_INCLUDE

// See loop3.hpp for example

namespace meta_math {

template <std::size_t Index0, std::size_t Max0, std::size_t Index1, std::size_t Max1>
struct loop2
{
    static const std::size_t index0= Index0 - 1, next_index0= Index0,
            	             index1= Index1 - 1, next_index1= Index1 + 1;
};


template <std::size_t Index0, std::size_t Max0, std::size_t Max1>
struct loop2<Index0, Max0, Max1, Max1>
{
    static const std::size_t index0= Index0 - 1, next_index0= Index0 + 1,
            	             index1= Max1 - 1, next_index1= 1;
};


template <std::size_t Max0, std::size_t Max1>
struct loop2<Max0, Max0, Max1, Max1>
{
    static const std::size_t index0= Max0 - 1,
            	             index1= Max1 - 1;
};


} // namespace meta_math

#endif // META_MATH_LOOP2_INCLUDE
