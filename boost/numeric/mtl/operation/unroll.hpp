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

#ifndef MTL_UNROLL_INCLUDE
#define MTL_UNROLL_INCLUDE

#include <boost/numeric/mtl/utility/enable_if.hpp>
#include <boost/numeric/mtl/vector/unrolled1.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>


namespace mtl {

/// Helper function for customizing loop unrolling in expressions
/** For instance,
    \code
      unroll<6>(u)= v + 3. * w; 
    \endcode
    will unroll the loop that computes
    \code
      u= v + 3. * w; 
    \endcode
    with block size 6.
**/
template <unsigned BSize, typename Coll> 
typename mtl::traits::enable_if_vector<Coll, vec::unrolled1<BSize, Coll> >::type
inline unroll(Coll& v)
{
    vampir_trace<7> tracer;
    return vec::unrolled1<BSize, Coll>(v);
}

} // namespace mtl

#endif // MTL_UNROLL_INCLUDE
