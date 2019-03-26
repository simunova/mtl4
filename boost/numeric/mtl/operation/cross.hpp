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

#ifndef MTL_VECTOR_CROSS_INCLUDE
#define MTL_VECTOR_CROSS_INCLUDE

#include <boost/numeric/mtl/concept/std_concept.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace mtl { namespace vec {

    namespace detail {

	// Result type of cross product
	template <typename Vector1, typename Vector2>
	struct cross_result
	{
	    typedef typename Multiplicable<typename Collection<Vector1>::value_type,
					   typename Collection<Vector2>::value_type>::result_type value;
	    typedef dense_vector<value>  type;
	};
    }

/// Cross product
/** Only exists for 3 and 7 dimensions.
    Consider specialization for fixed-size types
 **/
template <typename Vector1, typename Vector2>
typename detail::cross_result<Vector1, Vector2>::type
inline cross(const Vector1& v1, const Vector2& v2)
{
    vampir_trace<2002> tracer;
    MTL_THROW_IF((size(v1) != 3 && size(v1) != 7 ) || size(v1) != size(v2), incompatible_size());
    
    typename detail::cross_result<Vector1, Vector2>::type result(size(v1));

    if (size(v1) == 3) 
	for (unsigned i= 0; i < 3; i++) {
	    unsigned k= (i+1) % 3, l= (i+2) % 3;
	    result[i]= v1[k] * v2[l] - v1[l] * v2[k];
	}
    else  // must be 7 thus
	for (unsigned i= 0; i < 7; i++) {
	    unsigned k= (i+1) % 7, l= (i+3) % 7;
	    result[i]= v1[k] * v2[l] - v1[l] * v2[k];

	    k= (i+2) % 7, l= (i+6) % 7;
	    result[i]+= v1[k] * v2[l] - v1[l] * v2[k];

	    k= (i+4) % 7, l= (i+5) % 7;
	    result[i]+= v1[k] * v2[l] - v1[l] * v2[k];
	}    
    return result;
}


}} // namespace mtl::vector

#endif // MTL_VECTOR_CROSS_INCLUDE
