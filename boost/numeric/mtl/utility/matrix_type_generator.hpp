// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#ifndef MTL_MATRIX_TYPE_GENERATOR_INCLUDE
#define MTL_MATRIX_TYPE_GENERATOR_INCLUDE

#ifdef MTL_WITH_VARIADIC_TEMPLATE

#include <cstddef>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/map.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/void.hpp>

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/type_parameter.hpp>
#include <boost/numeric/mtl/utility/static_assert.hpp>
#include <boost/numeric/mtl/matrix/dimension.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>

namespace mtl {


    namespace type_para {

	template <typename DimPara>
	struct set_matrix_dimensions
	{
	    MTL_STATIC_ASSERT((! boost::is_same<DimPara, DimPara>::value), "Unsupported argument for matrix dimension.");
	};

	template <>
	struct set_matrix_dimensions<none>
	{
	    typedef mtl::non_fixed::dimensions type;
	};

	template <std::size_t ...Values>
	struct set_matrix_dimensions<dim<Values...> >
	{
	    MTL_STATIC_ASSERT((sizeof...(Values) == 2), "dim<rows, columns> must have exactly 2 arguments for matrices!");	    
	};

	template <std::size_t Rows, std::size_t Cols>
	struct set_matrix_dimensions<dim<Rows, Cols> >
	{
	    typedef mtl::fixed::dimensions<Rows, Cols> type;
	};

	
	template <typename TypePara>
	struct matrix_parameter_generator
	{
	    typedef typename boost::mpl::at<TypePara, orientation>::type            ori1;
	    typedef typename boost::mpl::if_<boost::is_same<ori1, none>,
					     row_major,
					     ori1>::type                            ori2;
	
	    typedef typename set_matrix_dimensions<typename boost::mpl::at<TypePara, dimensionality>::type>::type dim_type;
	    typedef typename boost::mpl::at<TypePara, size_type>::type                                            as_size;
	    typedef mtl::mat::parameters<
		ori2,
		index::c_index,
		dim_type,
		boost::is_same<typename boost::mpl::at<TypePara, location>::type, on_stack>::value,
		typename as_size::type
	    > type;
	};

	template <typename TypePara>
	struct matrix_default_density
	{
	    typedef typename boost::mpl::at<TypePara, density>::type init_density;
	    
	    typedef boost::mpl::map<
	    	boost::mpl::pair<banded, sparse>,       // Will be removed when dense banded is available
	    	boost::mpl::pair<compressed, sparse>,
	    	boost::mpl::pair<coordinate, sparse>,
	    	boost::mpl::pair<ellpack, sparse>,
	    	boost::mpl::pair<morton, dense>,
	    	boost::mpl::pair<none, dense>
	    >                                            default_map;

	    typedef typename boost::mpl::if_<
		boost::is_same<init_density, none>,
		typename boost::mpl::at<default_map, typename boost::mpl::at<TypePara, layout>::type>::type,
		init_density
	    >::type                                      type;
	};

	template <typename T>
	struct morton_matrix_mask
	{
	    MTL_STATIC_ASSERT(( !boost::is_same<T, T>::value), "Unknown type argument for Morton-order mask (internal error?).");
	};

	template <>
	struct morton_matrix_mask<none>
	{
	    static const std::size_t value= morton_mask;
	};

# ifndef _MSC_VER // creates problems on VS
	template <std::size_t ...Values>
	struct morton_matrix_mask<mask<Values...> >
	{
	    MTL_STATIC_ASSERT((sizeof...(Values) != 1), "Morton-order matrices must have exactly one mask");
	};	
# endif

	template <std::size_t Value>
	struct morton_matrix_mask<mask<Value> >
	{
	    static const std::size_t value= Value;
	};	


	// Generate dense matrix (Density should be dense or none)
	template <typename Value, typename Density, typename TypePara, typename MatrixPara>
	struct matrix_density_generator
	{
	    MTL_STATIC_ASSERT((boost::is_same<Density, dense>::value || boost::is_same<Density, none>::value), 
			      "Internal programm error.");
	    static const std::size_t my_mask= morton_matrix_mask<typename boost::mpl::at<TypePara, masking>::type>::value;
	    typedef boost::mpl::map<
	    	boost::mpl::pair<none, mat::dense2D<Value, MatrixPara> >,
		boost::mpl::pair<morton, mat::morton_dense<Value, my_mask, MatrixPara> >
	    >                                            type_map;
	    	    
	    typedef typename boost::mpl::at<type_map, typename boost::mpl::at<TypePara, layout>::type>::type type;
	    MTL_STATIC_ASSERT(( !boost::is_same<type, boost::mpl::void_>::value),
	    		       "The layout you providing cannot be used for dense matrices.");
	};


	// Generate sparse matrix
	template <typename Value, typename TypePara, typename MatrixPara>
	struct matrix_density_generator<Value, sparse, TypePara, MatrixPara>
	{
	    typedef boost::mpl::map<
	    	boost::mpl::pair<banded, mat::sparse_banded<Value, MatrixPara> >,
	    	boost::mpl::pair<compressed, mat::compressed2D<Value, MatrixPara> >,
	    	boost::mpl::pair<none, mat::compressed2D<Value, MatrixPara> >,
	    	boost::mpl::pair<coordinate, mat::coordinate2D<Value, MatrixPara> >,
	    	boost::mpl::pair<ellpack, mat::ell_matrix<Value, MatrixPara> >
	    >                                            type_map;
	    	    
	    typedef typename boost::mpl::at<type_map, typename boost::mpl::at<TypePara, layout>::type>::type type;
	    MTL_STATIC_ASSERT(( !boost::is_same<type, boost::mpl::void_>::value),
	    		       "The layout you providing cannot be used for sparse matrices.");
	};


	template <typename Value, typename TypePara>
	struct matrix_type_generator
	{
	    typedef typename matrix_parameter_generator<TypePara>::type  matrix_parameters;
	    // Imply density by layout
	    typedef typename matrix_default_density<TypePara>::type      my_density;
	    // Dispatch for density
	    typedef typename matrix_density_generator<Value, my_density, TypePara, matrix_parameters>::type type;
	};

    } // namespace type_para

# ifdef MTL_WITH_TEMPLATE_ALIAS
    /// Generated matrix type, see \ref type_generator
    template <typename Value, typename ...Parameters>
    using matrix= typename type_para::matrix_type_generator<Value, typename set_parameters<Parameters...>::type>::type;
# endif

} // namespace mtl

#endif // MTL_WITH_VARIADIC_TEMPLATE

#endif // MTL_MATRIX_TYPE_GENERATOR_INCLUDE
