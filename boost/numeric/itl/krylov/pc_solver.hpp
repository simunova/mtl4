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

#ifndef ITL_PC_SOLVER_INCLUDE
#define ITL_PC_SOLVER_INCLUDE

namespace itl {

#ifdef MTL_WITH_TEMPLATE_ALIAS

    template <typename Matrix, template <typename, typename, typename> class Solver,
	      template <typename> class Left, template <typename> class Right>
    using pc_solver= Solver<Matrix, Left<Matrix>, Right<Matrix> >;

#endif // MTL_WITH_TEMPLATE_ALIAS

} // namespace itl

#endif // ITL_PC_SOLVER_INCLUDE
