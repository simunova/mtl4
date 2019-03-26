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

#ifndef MTL_TRAITS_FLATCAT_INCLUDE
#define MTL_TRAITS_FLATCAT_INCLUDE

#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_base_of.hpp>

#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>

namespace mtl { namespace traits {


template <typename C, typename U1>
struct flatcat1_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    tag::universe
		   >::type
{};

template <typename T, typename U1>
struct flatcat1
  : flatcat1_c<typename category<T>::type, U1> {};

template <typename C, typename U1, typename U2>
struct flatcat2_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    flatcat1_c<C, U2>
		   >::type
{};

template <typename T, typename U1, typename U2>
struct flatcat2
  : flatcat2_c<typename category<T>::type, U1, U2> {};

template <typename C, typename U1, typename U2, typename U3>
struct flatcat3_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    flatcat2_c<C, U2, U3>
		   >::type
{};

template <typename T, typename U1, typename U2, typename U3>
struct flatcat3
  : flatcat3_c<typename category<T>::type, U1, U2, U3> {};


template <typename C, typename U1, typename U2, typename U3, typename U4>
struct flatcat4_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    flatcat3_c<C, U2, U3, U4>
		   >::type
{};

template <typename T, typename U1, typename U2, typename U3, typename U4>
struct flatcat4
  : flatcat4_c<typename category<T>::type, U1, U2, U3, U4> {};


template <typename C, typename U1, typename U2, typename U3, typename U4, typename U5>
struct flatcat5_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    flatcat4_c<C, U2, U3, U4, U5>
		   >::type
{};

template <typename T, typename U1, typename U2, typename U3, typename U4, typename U5>
struct flatcat5
  : flatcat5_c<typename category<T>::type, U1, U2, U3, U4, U5> {};

template <typename C, typename U1, typename U2, typename U3, typename U4, typename U5, typename U6>
struct flatcat6_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    flatcat5_c<C, U2, U3, U4, U5, U6>
		   >::type
{};

template <typename T, typename U1, typename U2, typename U3, typename U4, typename U5, typename U6>
struct flatcat6
  : flatcat6_c<typename category<T>::type, U1, U2, U3, U4, U5, U6> {};

template <typename C, typename U1, typename U2, typename U3, typename U4, typename U5, typename U6, typename U7>
struct flatcat7_c
  : boost::mpl::if_<boost::is_base_of<U1, C>,
		    tag::flat<U1>,
		    flatcat6_c<C, U2, U3, U4, U5, U6, U7>
		   >::type
{};

template <typename T, typename U1, typename U2, typename U3, typename U4, typename U5, typename U6, typename U7>
struct flatcat7
  : flatcat7_c<typename category<T>::type, U1, U2, U3, U4, U5, U6, U7> {};


// Some often reused flatcats:

#if 0
template <typename Matrix>
struct mat_cvec_flatcat
  : flatcat7<Matrix, tag::element_structure, tag::sparse_banded_matrix, tag::transposed_multi_vector, tag::hermitian_multi_vector, tag::multi_vector, tag::dense, tag::sparse> 
{};
#endif 

template <typename Matrix>
struct mat_cvec_flatcat
  : flatcat6<Matrix, tag::element_structure, tag::transposed_multi_vector, tag::hermitian_multi_vector, tag::multi_vector, tag::dense, tag::sparse> 
{};

template <typename Collection>
struct shape_flatcat
  : flatcat4<Collection, tag::matrix, tag::col_vector, tag::row_vector, tag::scalar>
{};

template <typename Collection>
struct sparsity_flatcat
  : flatcat2<Collection, tag::dense, tag::sparse>
{};

template <typename Collection>
struct cursor_flatcat
  : flatcat3<Collection, tag::has_fast_ra_cursor, tag::has_ra_cursor, tag::has_cursor>
{};

template <typename Collection>
struct iterator_flatcat
  : flatcat3<Collection, tag::has_fast_ra_iterator, tag::has_ra_iterator, tag::has_iterator>
{};

template <typename Collection>
struct has_iterator_flatcat
  : flatcat1<Collection, tag::has_iterator>
{};


template <typename Collection>
struct layout_flatcat
  : flatcat2<Collection, tag::has_2D_layout, tag::has_1D_layout>
{};


template <typename Matrix>
struct sub_matrix_flatcat
  : flatcat3<Matrix, tag::sub_divisible, tag::qsub_divisible, tag::has_sub_matrix>
{};

}} // namespace mtl::traits

#endif // MTL_TRAITS_FLATCAT_INCLUDE
