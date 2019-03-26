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

#ifndef ITL_PC_IC_0_INCLUDE
#define ITL_PC_IC_0_INCLUDE

#include <boost/mpl/bool.hpp>

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/operation/lower_trisolve.hpp>
#include <boost/numeric/mtl/operation/upper_trisolve.hpp>
#include <boost/numeric/mtl/matrix/upper.hpp>
#include <boost/numeric/mtl/matrix/strict_lower.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/itl/pc/solver.hpp>


namespace itl { namespace pc {

template <typename Matrix, typename Value= typename mtl::Collection<Matrix>::value_type>
class ic_0
{
  public:
    typedef Value                                                    value_type;
    typedef typename mtl::Collection<Matrix>::size_type              size_type;
    typedef ic_0                                                     self;

    typedef mtl::mat::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
    typedef mtl::mat::compressed2D<value_type, para>                      U_type;
#ifndef ITL_IC_0_ONE_MATRIX
    typedef U_type                                                   L_type;
#else
    typedef typename mtl::mat::transposed_view<U_type>            L_type;
#endif
    typedef mtl::mat::detail::lower_trisolve_t<L_type, mtl::tag::inverse_diagonal, true> lower_solver_t;
    typedef mtl::mat::detail::upper_trisolve_t<U_type, mtl::tag::inverse_diagonal, true> upper_solver_t;

    ic_0(const Matrix& A) : f(A, U), L(trans(U)), lower_solver(L), upper_solver(U) {}


    // solve x = U^* U y --> y= U^{-1} U^{-*} x
    template <typename Vector>
    Vector solve(const Vector& x) const
    {
	mtl::vampir_trace<5036> tracer;
	return inverse_upper_trisolve(U, inverse_lower_trisolve(adjoint(U), x));
    }

    // solve x = U^* y --> y0= U^{-*} x
    template <typename VectorIn, typename VectorOut>
    const VectorOut& solve_lower(const VectorIn& x, VectorOut&) const
    {
	static VectorOut y0;
	y0.change_resource(resource(x));
	lower_solver(x, y0);
	return y0;
    }

    // solve x = U^* U y --> y= U^{-1} U^{-*} x
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5037> tracer;
	const VectorOut& y0= solve_lower(x, y);

	y.checked_change_resource(x);
	upper_solver(y0, y);
    }

    // solve x = (LU)^* y --> y= L^{-*} U^{-*} x
    template <typename Vector>
    Vector adjoint_solve(const Vector& x) const
    {
	mtl::vampir_trace<5044> tracer;
	return solve(x);
    }

    // solve x = (LU)^* y --> y= L^{-*} U^{-*} x
    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5044> tracer;
	solve(x, y); 
    }


    L_type get_L() { return L_type(L); }
    U_type get_U() { return U; }

  protected:
    template <typename VectorOut, typename Solver> friend struct ic_0_evaluator;

    // Dummy type to perform factorization in initializer list not in 
    struct factorizer
    {
	factorizer(const Matrix &A, U_type& U)
	{   factorize(A, U, mtl::traits::is_sparse<Matrix>(), boost::is_same<Value, typename mtl::Collection<Matrix>::value_type>());  }

	template <typename T>
	void factorize(const Matrix&, U_type&, boost::mpl::false_, T)
	{   MTL_THROW_IF(true, mtl::logic_error("IC(0) is not suited for dense matrices"));	}

	// When we change the value_type then the factorization is still performed with that of A
	template <typename UF>
	void factorize(const Matrix& A, UF& U, boost::mpl::true_, boost::mpl::false_)
	{
	    typedef mtl::mat::compressed2D<typename mtl::Collection<Matrix>::value_type, para> tmp_type;
	    tmp_type U_tmp;
	    factorize(A, U_tmp, boost::mpl::true_(), boost::mpl::true_());
	    U= U_tmp;
	}

	// Factorization adapted from Saad
	// Undefined (runtime) behavior if matrix is not symmetric 
	// UF is type for the factorization
	template <typename UF> 
	void factorize(const Matrix& A, UF& U, boost::mpl::true_, boost::mpl::true_)
	{
	    using namespace mtl; using namespace mtl::tag;  using mtl::traits::range_generator;  
	    using math::reciprocal; using mtl::mat::upper;
	    mtl::vampir_trace<5035> tracer;

	    // For the factorization we take still the value_type of A and later we copy it maybe to another value_type
	    typedef typename mtl::Collection<Matrix>::value_type      value_type;
	    typedef typename range_generator<row, UF>::type       cur_type;    
	    typedef typename range_generator<nz, cur_type>::type      icur_type;            

	    MTL_THROW_IF(num_rows(A) != num_cols(A), mtl::matrix_not_square());
	    U= upper(A);

	    typename mtl::traits::col<UF>::type                   col(U);
	    typename mtl::traits::value<UF>::type                 value(U); 	

	    cur_type kc= begin<row>(U), kend= end<row>(U);
	    for (size_type k= 0; kc != kend; ++kc, ++k) {

		icur_type ic= begin<nz>(kc), iend= end<nz>(kc);
		MTL_DEBUG_THROW_IF(col(*ic) != k, mtl::missing_diagonal());

		// U[k][k]= 1.0 / sqrt(U[k][k]);
		value_type inv_dia= reciprocal(sqrt(value(*ic)));
		value(*ic, inv_dia);
		// icur_type jbegin= 
		++ic;
		for (; ic != iend; ++ic) {
		    // U[k][i] *= U[k][k]
		    value_type d= value(*ic) * inv_dia;
		    value(*ic, d);
		    size_type i= col(*ic);

		    // find non-zeros U[j][i] below U[k][i] for j in (k, i]
		    // 1. Go to ith row in U (== ith column in U)
		    cur_type irow(i, U); // = begin<row>(U); irow+= i;
		    // 2. Find nonzeros with col() in (k, i]
		    icur_type jc= begin<nz>(irow), jend= end<nz>(irow);
		    while (col(*jc) <= k)  ++jc;
		    while (col(*--jend) > i) ;
		    ++jend; 
		
		    for (; jc != jend; ++jc) {
			size_type j= col(*jc);
			U.lvalue(j, i)-= d * U[k][j];
		    }
		    // std::cout << "U after eliminating U[" << i << "][" << k << "] =\n" << U;
		}
	    }
	}
    };

    U_type                       U;
    factorizer                   f;
    L_type                       L;
    lower_solver_t               lower_solver;
    upper_solver_t               upper_solver;
}; 

#if 0
template <typename Matrix, typename Value, typename Vector>
struct ic_0_solver
  : mtl::assigner<ic_0_solver<Matrix, Value, Vector> >
{
    typedef ic_0<Matrix, Value> pc_type;

    ic_0_solver(const ic_0<Matrix, Value>& P, const Vector& x) : P(P), x(x) {}

    template <typename VectorOut>
    void assign_to(VectorOut& y) const
    {	P.solve(x, y);    }    

    const ic_0<Matrix, Value>& P; 
    const Vector&              x;
};
#endif

template <typename VectorOut, typename Solver>
struct ic_0_evaluator
{
    typedef typename Solver::pc_type                        pc_type;
    typedef typename pc_type::size_type                     size_type;
    typedef typename mtl::Collection<VectorOut>::value_type out_value_type;


    ic_0_evaluator(VectorOut& y, const Solver& s) 
      : y(y), s(s), U(s.P.U), y0(s.P.solve_lower(s.x, y)) { MTL_DEBUG_ARG(lr= 99999999); }


    void operator()(size_type i) { at<0>(i); }
    void operator[](size_type i) { at<0>(i); }

    template <unsigned Offset>
    void at(size_type r)
    {
#ifndef NDEBUG
	MTL_THROW_IF(r+Offset >= lr, mtl::logic_error("Traversal must be backward")); lr= r+Offset;
#endif
	size_type j0= U.ref_major()[r+Offset];
	const size_type cj1= U.ref_major()[r+Offset+1];
	MTL_DEBUG_THROW_IF(j0 == cj1 || U.ref_minor()[j0] != r+Offset, mtl::missing_diagonal());
	out_value_type rr= y0[r+Offset], dia= U.data[j0++];
	for (; j0 != cj1; ++j0) {
	    MTL_DEBUG_THROW_IF(U.ref_minor()[j0] <= r+Offset, mtl::logic_error("Matrix entries must be sorted for this."));
	    rr-= U.data[j0] * y[U.ref_minor()[j0]];
	}
	y[r+Offset]= rr * dia;
    }

    VectorOut&                               y;
    const Solver&                            s;
    const typename pc_type::U_type&          U;
    const VectorOut&                         y0;
    MTL_DEBUG_ARG(size_type                  lr;)
};

template <typename VectorOut, typename Solver>
inline std::size_t size(const ic_0_evaluator<VectorOut, Solver>& eval)
{   return size(eval.y); }

template <typename Matrix, typename Value, typename Vector>
solver<ic_0<Matrix, Value>, Vector, false>
inline solve(const ic_0<Matrix, Value>& P, const Vector& x)
{
    return solver<ic_0<Matrix, Value>, Vector, false>(P, x);
}

template <typename Matrix, typename Value, typename Vector>
solver<ic_0<Matrix, Value>, Vector, true>
inline adjoint_solve(const ic_0<Matrix, Value>& P, const Vector& x)
{
    return solver<ic_0<Matrix, Value>, Vector, true>(P, x);
}


}} // namespace itl::pc

namespace mtl { namespace vec {
    using itl::pc::size;
}} // namespace mtl::vector

#endif // ITL_PC_IC_0_INCLUDE
