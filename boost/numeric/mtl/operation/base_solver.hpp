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

#ifndef MTL_BASE_SOLVER_INCLUDE
#define MTL_BASE_SOLVER_INCLUDE

namespace mtl {

/// Helper class for common solver methods
class base_solver
{
  public:
    /// Default level is 1
    base_solver() : level(1) {}

    /// Set log level to \p level
    /** Concrete meaning of log level depends on specific solver.
	The general understanding is:
	- Level 0: no logging at all
	- Level 1: print final log (convergence and alike)
	- Level 2: print (some) intermediate log (e.g. final log of internal solver)
	- Level 3: more verbose (e.g. intermediate log of internal solver)
	- Level 4: even more verbose, ...   **/
    void set_log_level(unsigned level) { this->level= level; }

    /// Get log level
    unsigned log_level() const { return level; }

  private:
    unsigned level;

};

} // namespace mtl

#endif // MTL_BASE_SOLVER_INCLUDE
