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

#ifndef MTL_RANDOM_INCLUDE
#define MTL_RANDOM_INCLUDE

// Provisional implementation, do not use in production code!!!
// Will be reimplemented with boost::random

#include <cstdlib>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

template <typename T> 
struct seed 
{
    seed() { srand(17 * ++counter); }
    T operator()() const { return rand(); }

    static int counter;
};

template <typename T> 
int seed<T>::counter= 0;

namespace impl {

    // all non-distributed matrices except multi-vectors
    template <typename Coll, typename Generator>
    void inline random(Coll& A, Generator& generator, tag::universe, tag::matrix)
    {
	typedef typename Collection<Coll>::size_type size_type;
	mat::inserter<Coll> ins(A, A.dim2());
	for (size_type r= 0; r < num_rows(A); r++)
	    for (size_type c= 0; c < num_cols(A); c++)
		ins[r][c] << generator();
    }
    
    template <typename Coll, typename Generator>
    void inline random(Coll& A, Generator&, tag::universe, tag::multi_vector)
    {
	for (typename Collection<Coll>::size_type i= 0; i < num_cols(A); ++i)
	    random(A.vector(i));
    }

} // namespace impl


namespace vec {

    /// Fill vector with random values; generator must be a nullary function.
    template <typename Vector, typename Generator>
    typename mtl::traits::enable_if_vector<Vector>::type
    inline random(Vector& v, Generator& generator) 
    {
    vampir_trace<2021> tracer;
	typedef typename Collection<Vector>::size_type size_type;
	for (size_type i= 0; i < size(v); i++)
	    v[i]= generator();
    }

    /// Fill vector with random values.
    /** Currently done with rand(). Will be improved one day. You can provide
	your own generator as second argument. **/
    template <typename Vector>
    typename mtl::traits::enable_if_vector<Vector>::type
    inline random(Vector& v)
    {
	seed<typename Collection<Vector>::value_type> s;
	random(v, s);
    }


} // namespace vector

namespace mat {

    /// Fill matrix with random values; generator must be a nullary function.
    template <typename Matrix, typename Generator>
    typename mtl::traits::enable_if_matrix<Matrix>::type
    inline random(Matrix& A, Generator& generator) 
    {
    vampir_trace<4015> tracer;
	typename mtl::traits::category<Matrix>::type cat;
	mtl::impl::random(A, generator, cat, cat);
    }

    /// Fill matrix with random values.
    /** Currently done with rand(). Will be improved one day. You can provide
	your own generator as second argument. **/
    template <typename Matrix>
    typename mtl::traits::enable_if_matrix<Matrix>::type
    inline random(Matrix& A) 
    {
	seed<typename Collection<Matrix>::value_type> s; 
	random(A, s);
    }

} // namespace matrix

} // namespace mtl

#endif // MTL_RANDOM_INCLUDE
