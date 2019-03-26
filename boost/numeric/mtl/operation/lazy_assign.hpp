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

#ifndef MTL_LAZY_ASSIGN_INCLUDE
#define MTL_LAZY_ASSIGN_INCLUDE

namespace mtl {


/// Helper class for lazy assignment semantics
template <typename T, typename U, typename Assign>
struct lazy_assign
{
    typedef Assign  assign_type;

    lazy_assign(T& first, const U& second) : first(first), second(second) {} 
    void delay_assign() const {}    

    T&       first;
    const U& second;

};

} // namespace mtl

#endif // MTL_LAZY_ASSIGN_INCLUDE
