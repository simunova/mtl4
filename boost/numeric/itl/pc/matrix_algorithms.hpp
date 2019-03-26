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
//
// Algorithm inspired by Nick Vannieuwenhoven, written by Cornelius Steinhardt



#ifndef MTL_MATRIX_ALGORITHMS_INCLUDE
#define MTL_MATRIX_ALGORITHMS_INCLUDE

#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/matrix/element.hpp>
#include <boost/numeric/mtl/matrix/element_structure.hpp>

#include <iostream>

namespace mtl {  namespace mat {

/// Construct the sparse data structure from an elementstructure 
template< typename ElementStructure, typename Matrix, typename Vector> 
void assemble_compressed(const ElementStructure& es, Matrix& A, Vector& order) 
{

    typedef typename ElementStructure::element_type::value_type   value_type;
    typedef typename ElementStructure::element_iterator           iterator;
    typedef typename ElementStructure::element_type               element_type;
    typedef typename element_type::index_type                     index_type;
    typedef typename element_type::matrix_type                    matrix_type;
    typedef typename matrix_type::size_type                       size_type;
    A.change_dim(es.get_total_vars(), es.get_total_vars());
    set_to_zero(A);
    value_type zero(0);
	
    {//start inserterblock
	mtl::mat::inserter<Matrix, mtl::operations::update_plus<value_type> >  ins(A);
	for (iterator it = es.element_begin(); it != es.element_end(); ++it) {
	    element_type& element = *it;
	    const index_type& idx = element.get_indices();
	    matrix_type& values = element.get_values();
	    for (std::size_t i = 0; i < element.nb_vars(); ++i) 
		for (std::size_t j = 0; j < element.nb_vars(); ++j) 
		    if (values(i,j) != zero) 
			ins[size_type(order(idx(i)))][size_type(order(idx(j)))] << values(i,j);		    
	}
    }//end inserterblock
}


/// Construct the sparse data structure from an elementstructure 
template< typename ElementStructure, typename Matrix> 
void assemble_compressed(const ElementStructure& es, Matrix& A) 
{
  typedef typename  ElementStructure::element_type::matrix_type::size_type size_type;
  mtl::dense_vector<size_type> ident(es.get_total_vars());
  iota(ident);
  assemble_compressed(es, A, ident);
}

}}//end namespace mtl

#endif // MTL_MATRIX_ALGORITHMS_INCLUDE
