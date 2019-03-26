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

#ifndef MTL_DIAGONAL_INCLUDE
#define MTL_DIAGONAL_INCLUDE


#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

    namespace vec {

	/// Transform a vector into a diagonal matrix
	template <typename Vector>
	mtl::mat::compressed2D<typename Collection<Vector>::value_type, mat::parameters<> >
	// typename mtl::traits::enable_if_vector<Vector, mtl::mat::compressed2D<typename Collection<Vector>::value_type> >::type
	inline diagonal(const Vector& v)
	{
	    vampir_trace<2016> tracer;
	    typedef mtl::mat::compressed2D<typename Collection<Vector>::value_type, mat::parameters<> > matrix_type;
	    matrix_type                           D(size(v), size(v));
	    D= 0;
	    mtl::mat::inserter<matrix_type>    ins(D, 1);
	    for (typename Collection<Vector>::size_type i= 0; i < size(v); ++i)
		ins[i][i] << v[i];

	    return D;
	}
    } 

    namespace mat {

	/// Return the vector with the diagonal of the matrix
	template <typename Matrix>
	// typename mtl::traits::enable_if_matrix<Matrix, conj_view<Matrix> >::type
	mtl::vec::dense_vector<typename Collection<Matrix>::value_type, vec::parameters<> >
	inline diagonal(const Matrix& A)
	{
	    vampir_trace<3007> tracer;
	    using std::min;
	    typedef typename Collection<Matrix>::size_type size_type;
	    size_type n= min(num_rows(A), num_cols(A));
		mtl::vec::dense_vector<typename Collection<Matrix>::value_type, vec::parameters<> > v(n);

	    for (size_type i= 0; i < n; ++i)
		v[i]= A[i][i];
	    return v;
	}
    } 

} // namespace mtl

#endif // MTL_DIAGONAL_INCLUDE
