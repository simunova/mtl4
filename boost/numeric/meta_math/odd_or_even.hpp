// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Shikhar Vashistha
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef META_MATH_ODD_OR_EVEN
#define META_MATH_ODD_OR_EVEN

namespace meta_math{

template <long long n>
struct odd_or_even{
        static long long const value = (n^1)&1==true ? 0 : 1;
        // 0 corresponds to even number and 1 corresponds to odd number.
        };
    }// namespace meta_math

#endif// META_MATH_ODD_OR_EVEN

