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

#ifndef MTL_MATRIX_MULTI_VECTOR_RANGE_INCLUDE
#define MTL_MATRIX_MULTI_VECTOR_RANGE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/matrix/multi_vector.hpp>

namespace mtl { namespace mat {

// So far only const, might be refactored for mutability later
template <typename Vector>
class multi_vector_range
{
  public:
    typedef multi_vector_range                       self;
    typedef typename Collection<Vector>::size_type   size_type;
    typedef typename Collection<Vector>::value_type  value_type;
    typedef const value_type&                        const_reference;

    multi_vector_range(const multi_vector<Vector>& ref, const irange& r) : ref(ref), r(r) {}

    const_reference operator() (size_type i, size_type j) const { return ref[r.to_range(j)][i]; }
    const Vector& vector(size_type i) const { return ref.vector(r.to_range(i)); }

    /// Number of rows
    friend size_type num_rows(const self& A) { return num_rows(A.ref); }

    /// Number of columns
    friend size_type num_cols(const self& A) { return A.r.size(); }

    /// Size as defined by number of rows times columns
    friend size_type size(const self& A) { return num_rows(A) * num_cols(A); }

  protected:  
    const multi_vector<Vector>& ref;
    const irange                r;
};


}} // namespace mtl::matrix

#endif // MTL_MATRIX_MULTI_VECTOR_RANGE_INCLUDE
