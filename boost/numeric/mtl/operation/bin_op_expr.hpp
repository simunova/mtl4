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

#ifndef MTL_BIN_OP_EXPR_INCLUDE
#define MTL_BIN_OP_EXPR_INCLUDE

namespace mtl {

/// Minimalistic expression template for binary operation: keeps only references.
template <typename E1, typename E2>
struct bin_op_expr
{
    typedef bin_op_expr                          self;

    typedef E1                                   first_argument_type ;
    typedef E2                                   second_argument_type ;

    bin_op_expr( first_argument_type const& v1, second_argument_type const& v2 )
	: first( v1 ), second( v2 )
    {}

    first_argument_type const&  first ;
    second_argument_type const& second ;
};

} // namespace mtl

#endif // MTL_BIN_OP_EXPR_INCLUDE
