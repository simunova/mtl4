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

#ifndef MTL_TRAITS_OMP_SIZE_TYPE_INCLUDE
#define MTL_TRAITS_OMP_SIZE_TYPE_INCLUDE

#include <boost/mpl/if.hpp> 

namespace mtl { namespace traits {

# ifdef MTL_WITH_OPENMP

    template <typename T>
    struct omp_size_type
      : boost::mpl::if_c<(sizeof(T) > sizeof(int)), long int, int>
    {};

# else

    /// Type trait to provide size type w/wo OpenMP uniformely
    /** OpenMP emits warnings for unsigned ints and we therefor use signed ints only.
	Furthermore, we dispatch between int and long int. **/
    template <typename T>
    struct omp_size_type
    {
	typedef T type;
    };

# endif

}} // namespace mtl::traits

#endif // MTL_TRAITS_OMP_SIZE_TYPE_INCLUDE
