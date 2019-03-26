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

#ifndef ITL_PC_SOLVER_INCLUDE
#define ITL_PC_SOLVER_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/vector/assigner.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace itl { namespace pc {

/// Helper class for delayed (i.e. copy-free) evaluation of preconditioners
template <typename PC, typename Vector, bool adjoint= false>
struct solver
  : mtl::vec::assigner<solver<PC, Vector> >
{
    typedef PC  pc_type;

    /// Constructor taking preconditioner and source vector
    solver(const pc_type& P, const Vector& x) : P(P), x(x) {}

    /// Assign result to vector \p y, if possible without copying
    template <typename VectorOut>
    void assign_to(VectorOut& y) const
    {	
	mtl::vampir_trace<5055> tracer;
	assign_to(y, boost::mpl::bool_<adjoint>());    
    }    

  protected:

    template <typename VectorOut>
    void assign_to(VectorOut& y, boost::mpl::false_) const
    {	P.solve(x, y);    }    
    
    template <typename VectorOut>
    void assign_to(VectorOut& y, boost::mpl::true_) const
    {	P.adjoint_solve(x, y);    }    

    const pc_type&        P; 
    const Vector&         x;
};

}} // namespace itl::pc

#endif // ITL_PC_SOLVER_INCLUDE
