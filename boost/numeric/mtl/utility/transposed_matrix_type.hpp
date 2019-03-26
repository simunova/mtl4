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

#ifndef MTL_TRAITS_TRANSPOSED_MATRIX_TYPE_INCLUDE
#define MTL_TRAITS_TRANSPOSED_MATRIX_TYPE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/transposed_orientation.hpp>

namespace mtl { namespace traits {

template <class T> struct transposed_matrix_parameter {};

template <typename O, typename I, typename D, bool S, typename ST>
struct transposed_matrix_parameter<mat::parameters<O, I, D, S, ST> >
{
    typedef mat::parameters<typename transposed_orientation<O>::type, I, D, S, ST>  type;
};

template <class T> struct transposed_matrix_type {};

template <typename Value, typename Parameters>
struct transposed_matrix_type<mat::dense2D<Value, Parameters> >
{
    typedef mat::dense2D<Value, typename transposed_matrix_parameter<Parameters>::type> type;
};

template <typename Value, typename Parameters>
struct transposed_matrix_type<mat::compressed2D<Value, Parameters> >
{
    typedef mat::compressed2D<Value, typename transposed_matrix_parameter<Parameters>::type> type;
};

template <typename Value, std::size_t Mask, typename Parameters>
struct transposed_matrix_type<mat::morton_dense<Value, Mask, Parameters> >
{
    typedef mat::morton_dense<Value, Mask, typename transposed_matrix_parameter<Parameters>::type> type;
};




template <class T> struct transposed_sparse_matrix_type {};

template <typename Value, typename Parameters>
struct transposed_sparse_matrix_type<mat::compressed2D<Value, Parameters> >
{
    typedef mat::compressed2D<Value, typename transposed_matrix_parameter<Parameters>::type> type;
};

template <typename Matrix>
struct transposed_sparse_matrix_type<mat::banded_view<Matrix> >
{
    typedef typename transposed_sparse_matrix_type<Matrix>::type type;
};


template <typename Value, typename Parameters>
struct transposed_sparse_matrix_type<mat::transposed_view<mat::compressed2D<Value, Parameters> > >
{
    typedef mat::compressed2D<Value, Parameters> type;
};

template <typename Value, typename Parameters>
struct transposed_sparse_matrix_type<mat::transposed_view<const mat::compressed2D<Value, Parameters> > >
{
    typedef mat::compressed2D<Value, Parameters> type;
};

}} // namespace mtl::traits

#endif // MTL_TRAITS_TRANSPOSED_MATRIX_TYPE_INCLUDE
