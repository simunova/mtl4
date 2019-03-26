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

#ifndef MTL_PRINT_SIZE_INCLUDE
#define MTL_PRINT_SIZE_INCLUDE

#include <sstream>

#include <boost/numeric/mtl/operation/look_at_each_nonzero.hpp>

namespace mtl { 

    struct print_size_max
    {
	print_size_max() : max(0) {}

	template <typename T>
	void operator()(const T& x) 
	{ 
	    std::ostringstream st;
	    st << x;
	    std::size_t s= st.str().size();
	    // std::cout << "Size is " << s << ": " << st.str() << '\n';
	    if (s > max) max= s; 
	}

	std::size_t max;
    };

    namespace mat {
	    
	template <typename Matrix>
	std::size_t inline print_size(const Matrix& A)
	{
	    print_size_max p;
	    // look_at_each_nonzero(A, p); // Doesn't work on all expressions yet
	    for (size_t r = 0, nr = num_rows(A); r < nr; ++r)
		for (size_t c = 0, nc= num_cols(A); c < nc; ++c) 
		    p(A(r, c));
	    return p.max;
	}
    }

} // namespace mtl

#endif // MTL_PRINT_SIZE_INCLUDE
