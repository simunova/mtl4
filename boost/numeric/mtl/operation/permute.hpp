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

#ifndef MTL_PERMUTE_INCLUDE
#define MTL_PERMUTE_INCLUDE

#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>

namespace mtl {

namespace vec {

    /// Permute the entries of \p v according to permutation \p p.
    template <typename PermutationVector, typename Vector>
    typename mtl::traits::enable_if_vector<Vector, Vector>::type
    permute(const PermutationVector& p, const Vector& v)
    {
        std::size_t s = size(v);
        Vector result(s);
        for (std::size_t i = 0; i < s; ++i) {
            std::size_t j = p[i];
            if (j >= s)
                throw index_out_of_range("Index in permutation vector out of range (w.r.t. permuted vector).");
            result[j] = v[i];
        }
        return result;
    }

    /// Permute reversely the entries of \p v according to permutation \p p, inverse operation of permute.
    template <typename PermutationVector, typename Vector>
    typename mtl::traits::enable_if_vector<Vector, Vector>::type
    reverse_permute(const PermutationVector& p, const Vector& v)
    {
        std::size_t s = size(v);
        Vector result(s);
        for (std::size_t i = 0; i < s; ++i) {
            std::size_t j = p[i];
            if (j >= s)
                throw index_out_of_range("Index in permutation vector out of range (w.r.t. permuted vector).");
            result[i] = v[j];
        }
        return result;
    }


}                                                           // namespace vec

} // namespace mtl

#endif // MTL_PERMUTE_INCLUDE
