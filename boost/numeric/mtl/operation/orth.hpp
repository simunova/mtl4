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

#ifndef MTL_ORTH_INCLUDE
#define MTL_ORTH_INCLUDE

#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/operation/size1D.hpp>
#include <boost/numeric/mtl/operation/entry1D.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {

    namespace impl {

	template <typename VVector>
	inline void orth(VVector& v, typename mtl::Collection<VVector>::size_type j, tag::vector)
	{
		vampir_trace<2018> tracer;
	    using mtl::two_norm; using mtl::size1D;
	    MTL_DEBUG_THROW_IF(is_negative(j) || j >= size1D(v), index_out_of_range());

	    typedef typename mtl::Collection<VVector>::size_type  Size;
	    for (Size i= 0; i < j; ++i)
		entry1D(v, j)-= dot(entry1D(v, i), entry1D(v, j)) * entry1D(v, i);
	    entry1D(v, j)/= two_norm(entry1D(v, j));
	}

	template <typename VVector>
	inline void orth(VVector& v, tag::vector)
	{
	    typedef typename mtl::Collection<VVector>::size_type  Size;
	    using mtl::size1D;
	    for (Size j= 0; j < size1D(v); ++j)
		orth(v, j, tag::vector());
	}


	template <typename VVector>
	mtl::mat::dense2D<typename mtl::Collection
		   <typename mtl::Collection<VVector>::value_type
		    >::value_type, mat::parameters<> >
	inline orthogonalize_factors(VVector& v, tag::vector)
	{	
	    vampir_trace<2019> tracer;
	    using ::mtl::two_norm; using math::zero; using mtl::size1D;
	    typedef typename mtl::Collection<VVector>::size_type  Size;
	    typedef typename mtl::Collection<VVector>::value_type Vector;
	    typedef typename mtl::Collection<Vector>::value_type  Scalar;

	    mtl::mat::dense2D<Scalar, mat::parameters<> > tau(size1D(v), size1D(v));
	    tau= zero(Scalar());

	    if (size1D(v) == 0)
		return tau;

	    tau[0][0]= dot(entry1D(v, 0), entry1D(v, 0));

	    for (Size j= 1; j < size1D(v); ++j) {
#ifdef MTL_WITH_FUSED_ORTHOGONALIZATION
		Scalar t= dot(entry1D(v, 0), entry1D(v, j)) / tau[0][0], t2;
		tau[0][j]= t;
		for (Size i= 1; i < j; ++i) {
		    (lazy(entry1D(v, j))-= t * entry1D(v, i-1)) || (lazy(t2)= lazy_dot(entry1D(v, i), entry1D(v, j)));
		    t= tau[i][j]= t2 / tau[i][i];
		}
		entry1D(v, j)-= t * entry1D(v, j-1);
#else
		for (Size i= 0; i < j; ++i) {
		    Scalar t= dot(entry1D(v, i), entry1D(v, j)) / tau[i][i];
		    tau[i][j]= t;
		    entry1D(v, j)-= t * entry1D(v, i);
		}
#endif
		tau[j][j]= dot(entry1D(v, j), entry1D(v, j));
	    }
	    return tau;
	}

    } // impl



/*! Orthonormalize a vector of vectors.

    The outer type must be a random access collection and
    the vector type must provide a dot function. 
    For instance dense_vector<dense_vector<double> > or
    std::vector<dense_vector<std::complex<double> > > are eligible.
    It is planned to implement the function for matrices as well
    where the columns will be ortho-normalized.
**/
template <typename Value>
inline void orth(Value& value)
{
    impl::orth(value, typename traits::category<Value>::type());
}

/*! Orthonormalize the i-th entry of a vector of vectors.

    The i-th vector is orthogonalized w.r.t. to the preceeding ones and
    consecutively normalized.
    The outer type must be a random access collection and
    the vector type must provide a dot function. 
    For instance dense_vector<dense_vector<double> > or
    std::vector<dense_vector<std::complex<double> > > are eligible.
    It is planned to implement the function for matrices as well
    where the columns will be ortho-normalized.
**/
template <typename Value>
inline void orth(Value& value, typename mtl::Collection<Value>::size_type i)
{
    impl::orth(value, i, typename traits::category<Value>::type());
}


/*! Orthogonalize a vector of vectors.

    Opposed to orth the vectors are not normalized. 
    An upper matrix with the factors used in the orthogonalization is returned.
    The diagonal contains dot(v[i], v[i]).
    The returned factors are for instance used in bicgstab_ell.
    The outer type must be a random access collection and
    the vector type must provide a dot function. 
    For instance dense_vector<dense_vector<double> > or
    std::vector<dense_vector<std::complex<double> > > are eligible.
**/
template <typename Value>
mtl::mat::dense2D<typename mtl::Collection
	<typename mtl::Collection<Value>::value_type
	 >::value_type, mat::parameters<>  >
inline orthogonalize_factors(Value& v)
{
    return impl::orthogonalize_factors(v, typename traits::category<Value>::type());
}

} // namespace vector

namespace mat {

    // If other matrix types will be supported within a template function, it needs reimplementation!!!
    template <typename Vector>
    inline void orth(multi_vector<Vector>& A)
    {
	mtl::vec::impl::orth(A, mtl::tag::vector());
    }
}

using vec::orth;
using vec::orthogonalize_factors;

} // namespace mtl

#endif // MTL_ORTH_INCLUDE
