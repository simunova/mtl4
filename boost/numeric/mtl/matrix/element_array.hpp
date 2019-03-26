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

#ifndef MTL_ELEMENT_ARRAY_INCLUDE
#define MTL_ELEMENT_ARRAY_INCLUDE

namespace mtl { namespace mat {

template <typename Array, typename Rows, typename Cols>
struct element_array_t
{
    explicit element_array_t(const Array& array, const Rows& rows, const Cols& cols)
	: array(array), rows(rows), cols(cols)
    {}
    
    const Array& array;
    const Rows&  rows;
    const Cols&  cols;
};


template <typename Array, typename Rows, typename Cols>
element_array_t<Array, Rows, Cols>
inline element_array(const Array& array, const Rows& rows, const Cols& cols)
{
    return element_array_t<Array, Rows, Cols>(array, rows, cols);
}

template <typename Array, typename Rows>
element_array_t<Array, Rows, Rows>
inline element_array(const Array& array, const Rows& rows)
{
    return element_array_t<Array, Rows, Rows>(array, rows, rows);
}

} // namespace matrix

using mat::element_array;

} // namespace mtl

#endif // MTL_ELEMENT_ARRAY_INCLUDE
