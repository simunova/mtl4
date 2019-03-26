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

#ifndef MTL_SET_TO_0_INCLUDE
#define MTL_SET_TO_0_INCLUDE

#include <algorithm>
#include <cassert>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/ashape.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl {

    // Forward declarations
    namespace mat { 
	template <typename Coll> 
	typename mtl::traits::enable_if_matrix<Coll>::type 
	set_to_zero(Coll& collection); 
    }
    namespace vec { 
	template <typename Coll> 
	typename mtl::traits::enable_if_vector<Coll>::type
	set_to_zero(Coll& collection); 
    }

    namespace impl {

	template <typename Coll>
	void set_to_zero(Coll& collection, tag::vector_ref, ashape::scal)
	{
	    using math::zero;
	    typename Collection<Coll>::value_type  ref, my_zero(zero(ref));
	    for (typename Collection<Coll>::size_type i= 0; i < size(collection); ++i)
		collection[i]= my_zero;
	}

	template <typename Coll>
	void set_to_zero(Coll& collection, tag::contiguous_dense, ashape::scal)
	{
	    using math::zero;
	    typename Collection<Coll>::value_type  ref, my_zero(zero(ref));

	    std::fill(collection.elements(), collection.elements()+collection.used_memory(), my_zero);
	}

	template <typename Coll>
	void set_to_zero(Coll& collection, tag::std_vector, ashape::scal)
	{
	    using math::zero;
	    typename Collection<Coll>::value_type  ref, my_zero(zero(ref));

	    std::fill(collection.begin(), collection.end(), my_zero);
	}

	template <typename Matrix>
	void set_to_zero(Matrix& matrix, tag::morton_dense, ashape::scal)
	{
	    using math::zero;
	    typename Collection<Matrix>::value_type  ref, my_zero(zero(ref));
	    // maybe faster to do it straight
	    // if performance problems we'll take care of the holes
	    // std::cout << "set_to_zero: used_memory = " << matrix.used_memory() << "\n";
	    std::fill(matrix.elements(), matrix.elements() + matrix.used_memory(), my_zero);

#if 0
	    for (int i= 0; i < matrix.num_rows(); i++)
	      for (int j= 0; j < matrix.num_cols(); j++)
		matrix[i][j]= my_zero;
#endif
	}	

	// For nested collection, we must consider the dimensions of the elements
	// (Morton-order is included in contiguous_dense)
	template <typename Coll>
	void set_to_zero(Coll& collection, tag::contiguous_dense, ashape::nonscal)
	{
	    for (typename Collection<Coll>::size_type i= 0; i < collection.used_memory(); ++i)
		set_to_zero(collection.value_n(i));
	}


	// Is approbriate for all sparse matrices and vectors (including collections as value_type)
	template <typename Coll>
	void set_to_zero(Coll& collection, tag::sparse, ashape::universe)
	{
	    collection.make_empty();
	}
	
	// Special treatment for multi_vector
	template <typename Coll>
	void set_to_zero(Coll& collection, tag::multi_vector, ashape::universe)
	{
	    using mtl::vec::set_to_zero;
	    for (typename Collection<Coll>::size_type i= 0; i < num_cols(collection); ++i)
		set_to_zero(collection.vector(i));
	}	
	
	template <typename Coll>
	bool has_strided_data(const Coll&) 
	{ return false; }

	template <typename Value, typename Parameter>
	bool has_strided_data(const mat::dense2D<Value, Parameter>& A)
	{ return A.has_strided_data(); }

	
	template <typename Matrix>
	void naive_set_to_zero(Matrix& A, tag::matrix, tag::dense)
	{
	    using math::zero;
	    typename Collection<Matrix>::value_type  ref, my_zero(zero(ref));

	    for (unsigned i= 0; i < num_rows(A); i++)
		for (unsigned j= 0; j < num_cols(A); j++)
		    A[i][j]= my_zero;
	}

	template <typename Matrix>
	void naive_set_to_zero(Matrix&, tag::matrix, tag::sparse)
	{
	    assert(true); // must not be called
	}

    }


namespace mat {

    /// Sets all values of a collection to 0
    /// More spefically the defined multiplicative identity element
    template <typename Coll>
    typename mtl::traits::enable_if_matrix<Coll>::type
    set_to_zero(Coll& collection)
    {
	using mtl::traits::category;
	vampir_trace<3031> tracer;
	typedef typename Collection<Coll>::value_type value_type;
	if (mtl::impl::has_strided_data(collection))
	    mtl::impl::naive_set_to_zero(collection, typename category<Coll>::type(), typename category<Coll>::type());
	else
	    mtl::impl::set_to_zero(collection, typename category<Coll>::type(),typename ashape::ashape<value_type>::type()); // 2. ashape ???
    }   
}

namespace vec {

    /// Sets all values of a collection to 0
    /// More spefically the defined multiplicative identity element
    template <typename Coll>
    typename mtl::traits::enable_if_vector<Coll>::type
    set_to_zero(Coll& collection)
    {
	using mtl::traits::category;
	vampir_trace<2029> tracer;
	typedef typename Collection<Coll>::value_type value_type;
	mtl::impl::set_to_zero(collection, typename category<Coll>::type(),typename ashape::ashape<value_type>::type());
    }

}

} // namespace mtl

#endif // MTL_SET_TO_0_INCLUDE
