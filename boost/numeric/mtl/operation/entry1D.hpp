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

#ifndef MTL_ENTRY1D_INCLUDE
#define MTL_ENTRY1D_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/matrix/multi_vector.hpp>

namespace mtl {

    namespace vec {

	template <typename Vector>
	inline typename Collection<Vector>::value_type const&
	entry1D(const Vector& v, typename Collection<Vector>::size_type i)
	{
	    return v[i];
	}

	template <typename Vector>
	inline typename Collection<Vector>::value_type& 
	entry1D(Vector& v, typename Collection<Vector>::size_type i)
	{
	    return v[i];
	}
    }

    namespace mat {

	template <typename Vector>
	inline Vector& entry1D(multi_vector<Vector>& A, typename Collection<Vector>::size_type i)
	{
	    return A.vector(i);
	}	

	template <typename Vector>
	inline Vector const& entry1D(const multi_vector<Vector>& A, typename Collection<Vector>::size_type i)
	{
	    return A.vector(i);
	}	
    }

} // namespace mtl

#endif // MTL_ENTRY1D_INCLUDE
