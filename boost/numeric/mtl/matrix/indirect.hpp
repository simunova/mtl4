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

#ifndef MTL_MATRIX_INDIRECT_INCLUDE
#define MTL_MATRIX_INDIRECT_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/detail/range_generator.hpp>
#include <boost/numeric/mtl/utility/complexity.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/iset.hpp>
#include <boost/numeric/mtl/operation/is_negative.hpp>
#include <boost/numeric/mtl/matrix/mat_expr.hpp>

namespace mtl { namespace mat {

/// Class for indirect access to matrix
/** So far only constant access, mutable access for appropriate matrices might be added later. **/
template <typename Matrix>
struct indirect
  : mat_expr<indirect<Matrix> >,
    const_crtp_base_matrix<indirect<Matrix>, typename Collection<Matrix>::value_type, typename Collection<Matrix>::size_type>
{
    typedef indirect                                self;
    typedef Matrix                                  other;
    typedef typename Collection<Matrix>::value_type value_type;
    typedef typename Collection<Matrix>::size_type  size_type;
    // if implementation uses const_reference -> change collection accordingly


    /// Construct from constant matrix reference and isets for rows and columns
    indirect(const Matrix& ref, const iset& rows, const iset& cols) : ref(ref), rows(rows), cols(cols) {}

    friend size_type inline num_rows(const self& A) { return A.rows.size(); } ///< Number of rows
    friend size_type inline num_cols(const self& A) { return A.cols.size(); } ///< Number of colums
    size_type nnz() const { return num_rows(*this) * num_cols(*this); } ///< Number of non-zeros 

    size_type dim1() const { return rows.size(); } ///< Dimension 1 is equal to number of rows
    size_type dim2() const { return cols.size(); } ///< Dimension 2 is equal to number of columns

    value_type operator() (size_type r, size_type c) const
    {   return (ref(rows[r], cols[c]));    } ///< Read A[r][c]

  private:
    const Matrix& ref;
    iset          rows, cols;
};

template <typename Matrix> 
inline std::size_t size(const indirect<Matrix>& A)
{     return num_rows(A) * num_rows(A); }

}} // namespace mtl::matrix


// -- Range generators in utility/range_generator.hpp
// -- Property maps in utility/property_maps.hpp


#endif // MTL_MATRIX_INDIRECT_INCLUDE
