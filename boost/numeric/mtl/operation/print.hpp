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

#ifndef MTL_PRINT_INCLUDE
#define MTL_PRINT_INCLUDE

#include <iostream>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/operation/print_size.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/print_vector.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>


namespace mtl { 

namespace detail {

    template <typename Collection>
    struct with_format_t
    {
	explicit with_format_t(const Collection& collection, int width, int precision) 
	    : collection(collection), width(width), precision(precision)
	{}
	
	const Collection& collection;
	int width, precision;
    };

    template <typename Collection>
    inline std::ostream& operator<< (std::ostream& out, with_format_t<Collection> const& value) 
    {
	return print(value.collection, out, value.width, value.precision);
    }
    

} // namespace detail

namespace mat {

    template <typename Matrix>
    inline std::ostream& operator<< (std::ostream& out, const mat_expr<Matrix>& expr)
    {
	const Matrix& A= static_cast<const Matrix&>(expr);
	return print_matrix(A, out, print_size(A), 0);
    }

    template <typename Value>
    inline std::ostream&
    print(Value const& value, std::ostream& out= std::cout, int width= 3, int precision= 2)
    {
	return print_matrix(value, out, width, precision);
    }

    template <typename Collection>
    inline mtl::detail::with_format_t<Collection> with_format(const Collection& collection, int width= 3, int precision= 2)
    {
	return mtl::detail::with_format_t<Collection>(collection, width, precision);
    }  
} // namespace matrix

namespace vec {

    
    template <typename Vector>
    inline std::ostream& operator<< (std::ostream& out, const vec::vec_expr<Vector>& expr)
    {
	return print_vector(static_cast<const Vector&>(expr), out, 0, 0);
    }


    template <typename Value>
    inline std::ostream&
    print(Value const& value, std::ostream& out= std::cout, int width= 3, int precision= 2)
    {
	return print_vector(value, out, width, precision);
    }

    template <typename Collection>
    inline mtl::detail::with_format_t<Collection> with_format(const Collection& collection, int width= 3, int precision= 2)
    {
	return mtl::detail::with_format_t<Collection>(collection, width, precision);
    }

} // namespace vector


} // mtl


#endif // MTL_PRINT_INCLUDE
