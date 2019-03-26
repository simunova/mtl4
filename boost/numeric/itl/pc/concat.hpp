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

#ifndef ITL_PC_CONCAT_INCLUDE
#define ITL_PC_CONCAT_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/itl/pc/solver.hpp>

namespace itl { namespace pc {

/// Class for concatenating \tparam PC1 and \tparam PC2
template <typename PC1, typename PC2, typename Matrix, bool Store1= true, bool Store2= true>
class concat
{
    typedef typename boost::mpl::if_c<Store1, PC1, const PC1&>::type pc1_type;
    typedef typename boost::mpl::if_c<Store2, PC2, const PC2&>::type pc2_type;

  public:
    /// Construct both preconditioners from matrix \p A
    explicit concat(const Matrix& A) : A(A), pc1(A), pc2(A)
    {
	BOOST_STATIC_ASSERT((Store1 && Store2));
    }

    /// Both preconditioners are already constructed and passed as arguments
    /** If pc1 or pc2 is only constructed temporarily in the constructor call,
	the according Store argument must be true; otherwise the preconditioner
	will be a stale reference.
	Conversely, if the preconditioner is already build outside the constructor call,
	the according Store argument should be false for not storing the preconditioner twice. **/
    concat(const Matrix& A, const PC1& pc1, const PC2& pc2) : A(A), pc1(pc1), pc2(pc2) {}

  private:
    template <typename VectorOut>
    VectorOut& create_r(const VectorOut& y) const
    {
	static VectorOut  r(resource(y));
	return r;
    }

    template <typename VectorOut>
    VectorOut& create_d(const VectorOut& y) const
    {
	static VectorOut  d(resource(y));
	return d;
    }

  public:
    /// Concatenated preconditioning: pc2 is applied regularly and pc1 afterwards by defect correction
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5058> tracer;
	y.checked_change_resource(x);
	pc2.solve(x, y);

	VectorOut &r= create_r(y), &d= create_d(y);
	r= x;
	r-= A * y;

	pc1.solve(r, d);
	y+= d;
    }


    /// Concatenated preconditioning: adjoint of pc1 is applied regularly and pc2's adjoint afterwards by defect correction
    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5059> tracer;
	y.checked_change_resource(x);
	pc1.adjoint_solve(x, y);

	VectorOut &r= create_r(y), &d= create_d(y);
	r= x;
	r-= adjoint(A) * y;

	pc2.adjoint_solve(r, d);
	y+= d;
    }

   private:
    const Matrix& A;
    pc1_type 	  pc1;
    pc2_type   	  pc2;
};

template <typename PC1, typename PC2, typename Matrix, bool Store1, bool Store2, typename Vector>
solver<concat<PC1, PC2, Matrix, Store1, Store2>, Vector, false>
inline solve(const concat<PC1, PC2, Matrix, Store1, Store2>& P, const Vector& x)
{
    return solver<concat<PC1, PC2, Matrix, Store1, Store2>, Vector, false>(P, x);
}

template <typename PC1, typename PC2, typename Matrix, bool Store1, bool Store2, typename Vector>
solver<concat<PC1, PC2, Matrix, Store1, Store2>, Vector, true>
inline adjoint_solve(const concat<PC1, PC2, Matrix, Store1, Store2>& P, const Vector& x)
{
    return solver<concat<PC1, PC2, Matrix, Store1, Store2>, Vector, true>(P, x);
}


}} // namespace itl::pc

#endif // ITL_PC_CONCAT_INCLUDE
