// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_IO_WRITE_AST_INCLUDE
#define MTL_IO_WRITE_AST_INCLUDE

#include <fstream>
#include <boost/numeric/mtl/io/write_ast_dispatch.hpp>

namespace mtl { namespace io {

/// write the abstract syntax tree (AST) to file name \p fname
template <typename Expr>
void write_ast(const Expr& expr, const char* fname)
{
    std::ofstream f(fname);
    f << "digraph G {\n  ordering = out;\n  edge [arrowhead=none];\n\n";
    write_ast_dispatch(expr, "t", f);

    f << "}\n";
}


}} // namespace mtl::io

#endif // MTL_IO_WRITE_AST_INCLUDE
