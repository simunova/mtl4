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

#ifndef MTL_VECTOR_MAKE_SPARSE_INCLUDE
#define MTL_VECTOR_MAKE_SPARSE_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/size.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl { namespace vec {

// Commands in Matlab
// S = sparse(i,j,s,m,n,nzmax)
// S = sparse(i,j,s,m,n)
// S = sparse(i,j,s)
// S = sparse(m,n)

template <typename SizeVector, typename ValueVector>
struct make_sparse_trait
{
    typedef typename Collection<SizeVector>::value_type  size_type;
    typedef typename Collection<ValueVector>::value_type value_type;
    typedef mat::parameters<row_major, index::c_index, mtl::non_fixed::dimensions, false, size_type> paras;
    typedef mat::compressed2D<value_type, paras>   type;
};

/// Generates an \p m by \p n matrix from the vectors \p rows, \p cols, and \p values.
/** A sparse matrix is created (compressed2D). The value type is the same as the element type
    of the value vector and the size type the same as the entries of the vectors with the row indices.
    Zero entries in \p values are ignored. Entries with same coordinates are added.
    Same as <a href="http://www.mathworks.de/de/help/matlab/ref/sparse.html">Matlab's sparse</a> function besides that it is zero-indexed.
 **/
template <typename SizeVector1, typename SizeVector2, typename ValueVector>
inline typename make_sparse_trait<SizeVector1, ValueVector>::type
make_sparse(const SizeVector1& rows, const SizeVector2& cols, const ValueVector& values,
	    std::size_t m, std::size_t n)
{
    MTL_THROW_IF(size(rows) != size(cols), incompatible_size());
    MTL_THROW_IF(size(rows) != size(values), incompatible_size());

    typedef make_sparse_trait<SizeVector1, ValueVector> traits;
    typedef typename traits::type       matrix_type;
    typedef typename traits::value_type value_type;
    typedef typename traits::size_type  size_type;

    size_type               ms= size_type(m), ns= size_type(n);  // shouldn't be needed :-!
    matrix_type             A(ms, ns);
    {
        mat::inserter<matrix_type, update_plus<value_type> > ins(A, size_type(size(rows) / m + 1));

	for (std::size_t i= 0; i < size(rows); i++)
	    if (values[i] != value_type(0))
		ins[rows[i]][cols[i]] << values[i];
    }
    return A;
}

/// Generates an \p m by \p n matrix from the vectors \p rows, \p cols, and \p values.
/** A sparse matrix is created (compressed2D). The value type is the same as the element type
    of the value vector and the size type the same as the entries of the vectors with the row indices.
    Zero entries in \p values are ignored. Entries with same coordinates are added. Last parameter
    is ignored and can be omitted.
    Same as <a href="http://www.mathworks.de/de/help/matlab/ref/sparse.html">Matlab's sparse</a> function besides that it is zero-indexed.
 **/
template <typename SizeVector1, typename SizeVector2, typename ValueVector>
inline typename make_sparse_trait<SizeVector1, ValueVector>::type
make_sparse(const SizeVector1& rows, const SizeVector2& cols, const ValueVector& values,
	    std::size_t m, std::size_t n, std::size_t)
{
    return make_sparse(rows, cols, values, m, n);
}

/// Generates a matrix from the vectors \p rows, \p cols, and \p values.
/** A sparse matrix is created (compressed2D). 
    The number of rows/columns is one plus the maximum of the entries of \p rows and \p cols.
    The value type is the same as the element type
    of the value vector and the size type the same as the entries of the vectors with the row indices.
    Zero entries in \p values are ignored. Entries with same coordinates are added.
    Same as <a href="http://www.mathworks.de/de/help/matlab/ref/sparse.html">Matlab's sparse</a> function besides that it is zero-indexed.
 **/
template <typename SizeVector1, typename SizeVector2, typename ValueVector>
inline typename make_sparse_trait<SizeVector1, ValueVector>::type
make_sparse(const SizeVector1& rows, const SizeVector2& cols, const ValueVector& values)
{
    return make_sparse(rows, cols, values, max(rows)+1, max(cols)+1);
}

/// Generates an empty \p m by \p n matrix.
/** A sparse matrix is created (compressed2D<double>). 
    Same as <a href="http://www.mathworks.de/de/help/matlab/ref/sparse.html">Matlab's sparse</a> function besides that it is zero-indexed.
 **/
inline mat::compressed2D<double> make_sparse(std::size_t m, std::size_t n)
{
    return mat::compressed2D<double>(m, n);
}





} // namespace :vector

    using vec::make_sparse;

} // namespace mtl

#endif // MTL_VECTOR_MAKE_SPARSE_INCLUDE
