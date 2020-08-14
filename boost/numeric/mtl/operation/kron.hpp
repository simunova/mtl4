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

#ifndef MTL_MAT_KRON_INCLUDE
#define MTL_MAT_KRON_INCLUDE

#include <cstddef>

#include <boost/numeric/mtl/matrix/dense2D.hpp>

// #include <boost/numeric/mtl/operation/set_to_zero.hpp>

namespace mtl { namespace mat {

/// Kronecker product (prototype implementation)
// A full-fledged version will be in MTL5 (with expression templates and allowing for arbitrary mixed types)    
template <typename Value, typename Para>
dense2D<Value, Para> kron(const dense2D<Value, Para>& A, const dense2D<Value, Para>& B)
{
    using std::size_t;
    size_t nra= num_rows(A), nca= num_cols(A), nrb= num_rows(B), ncb= num_cols(B);
        
    dense2D<Value, Para> P(nra*nrb, nca*ncb);
    for (size_t ra= 0; ra < nra; ++ra)  
        for (size_t ca= 0; ca < nca; ++ca)
            for (size_t rb= 0; rb < nrb; ++rb)  
                for (size_t cb= 0; cb < ncb; ++cb)
                    P[ra*nrb+rb][ca*ncb+cb]= A[ra][ca] * B[rb][cb];
    return P;
}
    
    

}} // namespace mtl::mat

#endif // MTL_MAT_KRON_INCLUDE
