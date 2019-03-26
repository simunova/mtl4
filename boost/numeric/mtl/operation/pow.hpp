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

#ifndef MTL_POW_INCLUDE
#define MTL_POW_INCLUDE

#include <boost/numeric/mtl/vector/map_view.hpp>

namespace mtl {
    
    namespace vec {

        /// Raise Vector \a v to power \a exp
        template <typename Vector, typename Exponent>
        pow_by_view<Vector, Exponent> pow(const Vector& v, const Exponent& exp)
        {
            return pow_by_view<Vector, Exponent>(v, exp);
        }
        
    } // namespace vec

} // namespace mtl

#endif // MTL_POW_INCLUDE
