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

#include <iostream>
#include <boost/type_traits.hpp>

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/detail/index.hpp>
#include <boost/numeric/mtl/utility/maybe.hpp>
#include <boost/numeric/mtl/operation/raw_copy.hpp>
#include <boost/numeric/mtl/operation/update.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>


using namespace std;



template <typename Orientation, typename Indexing>
void test_compressed2D_insertion()
{
    typedef mtl::mat::parameters<Orientation, Indexing, mtl::fixed::dimensions<8, 6> >         parameters;
    typedef mtl::compressed2D<int, parameters>                                              matrix_type;
    matrix_type   matrix;  

    typedef mtl::compressed2D<int, mtl::mat::parameters<Orientation> > matrix2_type;  
    matrix2_type E;
    {
	mtl::mat::inserter<matrix2_type> empty(E);
    }

    {
	mtl::mat::compressed2D_inserter<int, parameters>  i0(matrix, 3);
	i0(2, 2) << 6; i0(7, 2) << 17; 
    }
    {   // Inserter that overwrites the old values
	mtl::mat::compressed2D_inserter<int, parameters>  i1(matrix, 3);

	i1(0, 3) << 31; i1(3, 3) << 33; i1(6, 0) << 34 << 35; i1(4, 4) << 36 << 37;
    }

    cout << "\n\n";
    print_matrix(matrix); 
    MTL_THROW_IF(matrix(0, 3) != 31, mtl::runtime_error("Error overwriting empty value"));
    MTL_THROW_IF(matrix(3, 3) != 33, mtl::runtime_error("Error overwriting existing value"));
    MTL_THROW_IF(matrix(6, 0) != 35, mtl::runtime_error("Error overwriting empty value twice"));
    MTL_THROW_IF(matrix(4, 4) != 37, mtl::runtime_error("Error overwriting existing value twice"));

    {   // Inserter that adds to the old values
        mtl::mat::compressed2D_inserter<int, parameters, mtl::operations::update_plus<int> > i2(matrix, 3);    
 
	i2(2, 2) << 21; i2(2, 4) << 22; i2(6, 1) << 23; 
	i2(7, 2) << 24 << 2; i2(4, 2) << 25; i2(2, 5) << 26; 
	i2(0, 2) << 27; i2(3, 1) << 28; i2(4, 2) << 29; 
    }
    cout << "\n\n";
    print_matrix(matrix); 
    MTL_THROW_IF(matrix(0, 2) != 27, mtl::runtime_error("Error adding to empty value"));
    MTL_THROW_IF(matrix(2, 2) != 27, mtl::runtime_error("Error adding to existing value"));
    MTL_THROW_IF(matrix(4, 2) != 54, mtl::runtime_error("Error adding to existing value twice (in 2 statements)"));
    MTL_THROW_IF(matrix(7, 2) != 43, mtl::runtime_error("Error adding to existing value twice (in 1 statement)"));
    cout << "\n\n";
    {
	mtl::mat::inserter<matrix_type, mtl::operations::update_plus<int> >  i3(matrix, 7);
	i3(2, 2) << 1;
    }

    MTL_THROW_IF(matrix(2, 2) != 28, mtl::runtime_error("Error adding to existing value"));
}
 
int main(int, char**)
{
    test_compressed2D_insertion<mtl::row_major, mtl::index::c_index>();
    test_compressed2D_insertion<mtl::col_major, mtl::index::c_index>();

    return 0;
}
