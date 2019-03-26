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

#ifndef MTL_LOWER_TRISOLVE_INCLUDE
#define MTL_LOWER_TRISOLVE_INCLUDE

#include <boost/mpl/int.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/property_map.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/linear_algebra/inverse.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace mat {


namespace detail {

    /// Class that implements lower trisolver
    /** DiaTag can be tag::regular_diagonal, tag::unit_diagonal, or tag::inverse_diagonal.
	CompactStorage means that matrix contains only lower entries (strict lower when DiaTag == unit_diagonal). \sa \ref trisolve_object **/
    template <typename Matrix, typename DiaTag, bool CompactStorage= false>
    struct lower_trisolve_t
    {
	MTL_STATIC_ASSERT((boost::is_same<DiaTag, tag::regular_diagonal>::value
			   || boost::is_same<DiaTag, tag::unit_diagonal>::value
			   || boost::is_same<DiaTag, tag::inverse_diagonal>::value),
			  "DiaTag must be either tag::regular_diagonal, tag::unit_diagonal, or tag::inverse_diagonal.");

	typedef typename Collection<Matrix>::value_type         	       	    value_type;
	typedef typename Collection<Matrix>::size_type          	       	    size_type;
	typedef typename OrientedCollection<Matrix>::orientation	       	    my_orientation;
	typedef typename mtl::traits::category<Matrix>::type    	       	    my_category;
	typedef typename mtl::traits::range_generator<tag::major, Matrix>::type     a_cur_type; // row or col depending on Matrix    
	typedef typename mtl::traits::range_generator<tag::nz, a_cur_type>::type    a_icur_type;   

	/// Construction from matrix \p A
	lower_trisolve_t(const Matrix& A) : A(A), value_a(A), col_a(A), row_a(A)
	{    MTL_THROW_IF(num_rows(A) != num_cols(A), matrix_not_square());	}
	
	template <typename M, typename D, bool C>
	struct generic_version
	  : boost::mpl::int_<(mtl::traits::is_row_major<M>::value ? 0 : 2)
	                     + (boost::is_same<D, tag::unit_diagonal>::value ? 1 : 2)
                            > {};

	template <typename M, typename D, bool C>
	struct version
	  : generic_version<M, D, C> {};


	template <typename Value, typename Para, typename D>
	struct version<compressed2D<Value, Para>, D, true>
	  : boost::mpl::if_<mtl::traits::is_row_major<Para>,
			    typename boost::mpl::if_<boost::is_same<D, tag::unit_diagonal>,
						     boost::mpl::int_<5>,
						     boost::mpl::int_<6> >::type,
			    generic_version<compressed2D<Value, Para>, D, true>
	                   >::type {};

	/// Solve \p w = A * \p v
	template <typename VectorIn, typename VectorOut>
	void operator()(const VectorIn& v, VectorOut& w) const
	{   vampir_trace<5022> tracer; apply(v, w, version<Matrix, DiaTag, CompactStorage>()); }
	

      private:
	template <typename Value>
	Value inline lower_trisolve_diavalue(const Value& v, tag::regular_diagonal) const
	{   using math::reciprocal; return reciprocal(v); }

	template <typename Value>
	Value lower_trisolve_diavalue(const Value& v, tag::inverse_diagonal) const
	{  return v;	}    

	template <typename Tag> int dia_inc(Tag) { return 0; }
	int dia_inc(tag::unit_diagonal) { return 1; }

	// Generic row-major unit_diagonal
	template <typename VectorIn, typename VectorOut>
	void apply(const VectorIn& v, VectorOut& w, boost::mpl::int_<1>) const
	{
	    using namespace tag; 
	    a_cur_type ac= begin<row>(A), aend= end<row>(A); 
	    for (size_type r= 0; ac != aend; ++r, ++ac) {
		a_icur_type aic= begin<nz>(ac), aiend= CompactStorage ? end<nz>(ac) : lower_bound<nz>(ac, r);
		typename Collection<VectorOut>::value_type rr= v[r];
		for (; aic != aiend; ++aic) {
		    MTL_DEBUG_THROW_IF(col_a(*aic) >= r, logic_error("Matrix entries must be sorted for this."));
		    rr-= value_a(*aic) * w[col_a(*aic)];
		}
		w[r]= rr;
	    }
	}

	// Generic row-major not unit_diagonal
	template <typename VectorIn, typename VectorOut>
	void apply(const VectorIn& v, VectorOut& w, boost::mpl::int_<2>) const
	{
	    using namespace tag; 
	    a_cur_type ac= begin<row>(A), aend= end<row>(A); 
	    for (size_type r= 0; ac != aend; ++r, ++ac) {
		a_icur_type aic= begin<nz>(ac), aiend= CompactStorage ? end<nz>(ac) : lower_bound<nz>(ac, r+1);
		MTL_THROW_IF(aic == aiend, missing_diagonal());
		--aiend;
		MTL_THROW_IF(col_a(*aiend) != r, missing_diagonal());

		value_type dia= value_a(*aiend);
		typename Collection<VectorOut>::value_type rr= v[r];

		for (; aic != aiend; ++aic) {
		    MTL_DEBUG_THROW_IF(col_a(*aic) >= r, logic_error("Matrix entries must be sorted for this."));
		    rr-= value_a(*aic) * w[col_a(*aic)];
		}
		w[r]= rr * lower_trisolve_diavalue(dia, DiaTag());
	    }
	}	

	// Generic column-major unit_diagonal
	template <typename VectorIn, typename VectorOut>
	void apply(const VectorIn& v, VectorOut& w, boost::mpl::int_<3>) const
	{
	    using namespace tag; 
	    w= v;
	    a_cur_type ac= begin<col>(A), aend= end<col>(A); 
	    for (size_type r= 0; ac != aend; ++r, ++ac) {
		a_icur_type aic= CompactStorage ? begin<nz>(ac) : lower_bound<nz>(ac, r+1), aiend= end<nz>(ac);
		typename Collection<VectorOut>::value_type rr= w[r];

		for (; aic != aiend; ++aic) {
		    MTL_DEBUG_THROW_IF(row_a(*aic) <= r, logic_error("Matrix entries must be sorted for this."));
		    w[row_a(*aic)]-= value_a(*aic) * rr;
		}
	    }
	}

	// Generic column-major not unit_diagonal
	template <typename VectorIn, typename VectorOut>
	void apply(const VectorIn& v, VectorOut& w, boost::mpl::int_<4>) const
	{
	    using namespace tag;
	    w= v;
	    a_cur_type ac= begin<col>(A), aend= end<col>(A); 
	    for (size_type r= 0; ac != aend; ++r, ++ac) {
		a_icur_type aic= CompactStorage ? begin<nz>(ac) : lower_bound<nz>(ac, r), aiend= end<nz>(ac);
		MTL_DEBUG_THROW_IF(aic == aiend || row_a(*aic) != r, missing_diagonal());
		typename Collection<VectorOut>::value_type rr= w[r]*= lower_trisolve_diavalue(value_a(*aic), DiaTag());

		for (++aic; aic != aiend; ++aic) {
		    MTL_DEBUG_THROW_IF(row_a(*aic) <= r, logic_error("Matrix entries must be sorted for this."));
		    w[row_a(*aic)]-= value_a(*aic) * rr;
		}
	    }
	}

	// Tuning for IC_0 and similar using compressed2D row-major compact with implicit unit diagonal
	template <typename VectorIn, typename VectorOut>
	void apply(const VectorIn& v, VectorOut& w, boost::mpl::int_<5>) const
	{
	    vampir_trace<5048> tracer;
	    if (num_rows(A) == 0) return;
	    size_type j1= A.ref_major()[1];
	    for (size_type r= 0, rend= num_rows(A); r != rend; ++r) {
		size_type j0= j1; 
		j1= A.ref_major()[r+1];
		typename Collection<VectorOut>::value_type rr= v[r];
		for (; j0 != j1; ++j0) {
		    MTL_DEBUG_THROW_IF(A.ref_minor()[j0] > r, logic_error("Matrix entries from U in lower triangular."));
		    rr-= A.data[j0] * w[A.ref_minor()[j0]];
		}
		w[r]= rr;
	    }
	}	

	// Tuning for IC_0 and similar using compressed2D row-major compact with explicitly stored diagonal (possibly already inverted)
	template <typename VectorIn, typename VectorOut>
	void apply(const VectorIn& v, VectorOut& w, boost::mpl::int_<6>) const
	{
	    vampir_trace<5047> tracer;
	    for (size_type r= 0, rend= num_rows(A); r != rend; ++r) {
		size_type j0= A.ref_major()[r], j1= A.ref_major()[r+1];
		MTL_THROW_IF(j0 == j1, missing_diagonal());
		--j1;
		MTL_THROW_IF(A.ref_minor()[j1] != r, missing_diagonal());
		value_type dia= A.data[j1];
		typename Collection<VectorOut>::value_type rr= v[r];
		for (; j0 != j1; ++j0) {
		    MTL_DEBUG_THROW_IF(A.ref_minor()[j0] > r, logic_error("Matrix entries from U in lower triangular."));
		    rr-= A.data[j0] * w[A.ref_minor()[j0]];
		}
		w[r]= rr * lower_trisolve_diavalue(dia, DiaTag());
	    }
	}	

	const Matrix&                                    A;
	typename mtl::traits::const_value<Matrix>::type  value_a; 
	typename mtl::traits::col<Matrix>::type          col_a; 
	typename mtl::traits::row<Matrix>::type          row_a;
    };

}  // detail


/// Solves the lower triangular matrix A  with the rhs v and returns the solution vector
template <typename Matrix, typename Vector>
Vector inline lower_trisolve(const Matrix& A, const Vector& v)
{
    Vector w(resource(v));
    detail::lower_trisolve_t<Matrix, tag::regular_diagonal> solver(A); 
    solver(v, w);
    return w;
}

/// Solves the lower triangular matrix A  with the rhs v 
template <typename Matrix, typename VectorIn, typename VectorOut>
inline void lower_trisolve(const Matrix& A, const VectorIn& v, VectorOut& w)
{
    detail::lower_trisolve_t<Matrix, tag::regular_diagonal> solver(A); 
    solver(v, w);
}

/// Solves the lower triangular matrix A (only one's in the diagonal) with the rhs v and returns the solution vector
template <typename Matrix, typename Vector>
Vector inline unit_lower_trisolve(const Matrix& A, const Vector& v)
{
    Vector w(resource(v));
    detail::lower_trisolve_t<Matrix, tag::unit_diagonal> solver(A); 
    solver(v, w);
    return w;
}

/// Solves the lower triangular matrix A (only one's in the diagonal) with the rhs v and returns the solution vector
template <typename Matrix, typename VectorIn, typename VectorOut>
inline void unit_lower_trisolve(const Matrix& A, const VectorIn& v, VectorOut& w)
{
    detail::lower_trisolve_t<Matrix, tag::unit_diagonal> solver(A); 
    solver(v, w);
}

/// Solves the lower triangular matrix A (inverse the diagonal) with the rhs v and returns the solution vector
template <typename Matrix, typename Vector>
Vector inline inverse_lower_trisolve(const Matrix& A, const Vector& v)
{
    Vector w(resource(v));
    detail::lower_trisolve_t<Matrix, tag::inverse_diagonal> solver(A); 
    solver(v, w);
    return w;
}

/// Solves the lower triangular matrix A (inverse the diagonal) with the rhs v and returns the solution vector
template <typename Matrix, typename VectorIn, typename VectorOut>
inline void inverse_lower_trisolve(const Matrix& A, const VectorIn& v, VectorOut& w)
{
    detail::lower_trisolve_t<Matrix, tag::inverse_diagonal> solver(A); 
    solver(v, w);
}

template <typename Matrix, typename Vector, typename DiaTag>
Vector inline lower_trisolve(const Matrix& A, const Vector& v, DiaTag)
{
    Vector w(resource(v));
    detail::lower_trisolve_t<Matrix, DiaTag> solver(A); 
    solver(v, w);
    return w;
}

template <typename Matrix, typename VectorIn, typename VectorOut, typename DiaTag>
inline void lower_trisolve(const Matrix& A, const VectorIn& v, VectorOut& w, DiaTag)
{
    detail::lower_trisolve_t<Matrix, DiaTag> solver(A); 
    solver(v, w);
}

}} // namespace mtl::matrix

#endif // MTL_LOWER_TRISOLVE_INCLUDE
