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

#ifndef MTL_ENTRY_SIMILAR_INCLUDE
#define MTL_ENTRY_SIMILAR_INCLUDE



#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>


namespace mtl {

    namespace mat {
	
	template <typename Matrix, typename Value>
	bool inline entry_similar(const Matrix& A, typename Collection<Matrix>::size_type i, 
				  typename Collection<Matrix>::size_type j, const Value& v, const Value& eps, tag::universe)
	{
	    using std::abs;
	    return abs(A[i][j] - v) <= eps;
	}	

	/// Compares A[i][j] with \p v returns false if larger than eps
	/** Function works on distributed matrices where only one process tests and the other return true.
	    The functions enables writing tests that work in parallel (and all other platforms). **/
	template <typename Matrix, typename Value>
        bool inline entry_similar(const Matrix& A, typename Collection<Matrix>::size_type i, 
				  typename Collection<Matrix>::size_type j, const Value& v, const Value& eps)
	{
	    return entry_similar(A, i, j, v, eps,  typename mtl::traits::category<Matrix>::type());
	}

    } // namespacs matrix

    namespace vec {

	template <typename Vector, typename Value>
	bool inline entry_similar(const Vector& x, typename Collection<Vector>::size_type i, 
				  const Value& v, const Value& eps, tag::universe)
	{
	    using std::abs;
	    return abs(x[i] - v) <= eps;
	}	

	/// Compares x[i] with \p v returns false if larger than eps
	/** Function works on distributed matrices where only one process tests and the other return true.
	    The functions enables writing tests that work in parallel (and all other platforms). **/
	template <typename Vector, typename Value>
        bool inline entry_similar(const Vector& x, typename Collection<Vector>::size_type i, 
				  const Value& v, const Value& eps)
	{
	    return entry_similar(x, i, v, eps,  typename mtl::traits::category<Vector>::type());
	}

    }

} // namespace mtl

#endif // MTL_ENTRY_SIMILAR_INCLUDE
