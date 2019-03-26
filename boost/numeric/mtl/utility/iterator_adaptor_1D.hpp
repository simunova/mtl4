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

#ifndef MTL_ITERATOR_ADAPTOR_1D_INCLUDE
#define MTL_ITERATOR_ADAPTOR_1D_INCLUDE

#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/utility/iterator_adaptor.hpp>

namespace mtl { namespace traits {

// Range generator for one-dimensional iterators
// Implemented by means of cursors and value map

template <typename Matrix>
struct range_generator<tag::const_iter::nz, Matrix>
{
    typedef typename traits::const_value<Matrix>::type     map_type;
    typedef typename glas::tag::nz                      cursor_tag;
    typedef range_generator<cursor_tag, Matrix>            cursor_range;
    typedef typename cursor_range::type                    cursor_type;
    typedef typename Matrix::value_type                    value_type;
    typedef typename cursor_range::complexity              complexity;
    static int const                                       level = cursor_range::level;

    // iterator is type of range generator
    typedef utilities::const_iterator_adaptor<map_type, cursor_type, value_type>  type;

    type begin(const Matrix& matrix) const
    {
	return type(map_type(matrix), cursor_range().begin(matrix));
    }

    type end(const Matrix& matrix) const
    {
	return type(map_type(matrix), cursor_range().end(matrix));
    }
};


#if 0
template <typename Matrix>
struct range_generator<tag::iter::nz, Matrix>
{
    typedef typename traits::value<Matrix>::type                                  map_type;
    typedef typename glas::tag::nz                                             cursor_tag;
    typedef typename range_generator<cursor_tag, Matrix>::type                    cursor_type;
    typedef typename Matrix::value_type                                           value_type;
    typedef typename range_generator<cursor_tag, Matrix>::complexity              complexity;
    static int const level = range_generator<cursor_tag, Matrix>::level;

    // iterator is type of range generator
    typedef utilities::const_iterator_adaptor<map_type, cursor_type, value_type>  type;

    type begin(const Matrix& matrix) const
    {
	map_type map(matrix);
	return type(map, begin<cursor_tag>(matrix));
    }

    type end(const Matrix& matrix) const
    {
	map_type map(matrix);
	return type(map, end<cursor_tag>(matrix));
    }
};


template <typename Matrix>
struct range_generator<tag::const_iter::all, Matrix>
{
    typedef typename traits::const_value<Matrix>::type                            map_type;
    typedef typename glas::tag::all                                            cursor_tag;
    typedef typename range_generator<cursor_tag, Matrix>::type                    cursor_type;
    typedef typename Matrix::value_type                                           value_type;
    typedef typename range_generator<cursor_tag, Matrix>::complexity              complexity;
    static int const level = range_generator<cursor_tag, Matrix>::level;

    // iterator is type of range generator
    typedef utilities::const_iterator_adaptor<map_type, cursor_type, value_type>  type;

    type begin(const Matrix& matrix) const
    {
	return type(map_type(matrix), begin<cursor_tag>(matrix));
    }

    type end(const Matrix& matrix) const
    {
	return type(map_type(matrix), end<cursor_tag>(matrix));
    }
};

template <typename Matrix>
struct range_generator<tag::iter::all, Matrix>
{
    typedef typename traits::value<Matrix>::type                                  map_type;
    typedef typename glas::tag::all                                            cursor_tag;
    typedef typename range_generator<cursor_tag, Matrix>::type                    cursor_type;
    typedef typename Matrix::value_type                                           value_type;
    typedef typename range_generator<cursor_tag, Matrix>::complexity              complexity;
    static int const level = range_generator<cursor_tag, Matrix>::level;

    // iterator is type of range generator
    typedef utilities::const_iterator_adaptor<map_type, cursor_type, value_type>  type;

    type begin(const Matrix& matrix) const
    {
	map_type map(matrix);
	return type(map, begin<cursor_tag>(matrix));
    }

    type end(const Matrix& matrix) const
    {
	map_type map(matrix);
	return type(map, end<cursor_tag>(matrix));
    }
};


#endif






#if 0
template <typename Matrix>
struct morton_dense_row_const_iterator
    : utilities::const_iterator_adaptor<typename traits::const_value<Matrix>::type, morton_dense_row_cursor<Matrix::mask>,
					typename Matrix::value_type>
{
    static const unsigned long                          mask= Matrix::mask;
    typedef morton_dense_row_cursor<mask>               cursor_type;
    typedef typename traits::const_value<Matrix>::type  map_type;
    typedef typename Matrix::value_type                 value_type;
    typedef typename Matrix::size_type                  size_type;
    typedef utilities::iterator_adaptor<map_type, cursor_type, value_type> base;
    
    morton_dense_row_const_iterator(const Matrix& matrix, size_type row, size_type col)
	: base(map_type(matrix), cursor_type(row, col))
    {}
};

utilities::const_iterator_adaptor<typename traits::const_value<Matrix>::type, morton_dense_row_cursor<Matrix::mask>,
					typename Matrix::value_type>
#endif

}} // namespace mtl::traits

#endif // MTL_ITERATOR_ADAPTOR_1D_INCLUDE
