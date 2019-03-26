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

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/compressed2D.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/maybe.hpp>
#include <boost/numeric/mtl/utility/complexity.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>


using namespace std;

struct test_exception {};

template <typename Matrix>
void one_d_iteration(char const* name, Matrix & matrix)
{
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type                                 row(matrix);
    typename traits::col<Matrix>::type                                 col(matrix);
    typename traits::value<Matrix>::type                               value(matrix);
    typedef  mtl::tag::nz                                              tag;
    typedef typename traits::range_generator<tag, Matrix>::type        cursor_type;
    typedef typename traits::range_generator<tag, Matrix>::complexity  complexity;
    
    cout << name << "\nElements: " << complexity() << '\n';
    for (cursor_type cursor(mtl::begin<tag>(matrix)), cend(mtl::end<tag>(matrix)); cursor != cend; ++cursor) {
	cout << "matrix[" << row(*cursor) << ", " << col(*cursor) << "] = " << value(*cursor) << '\n';
	if (row(*cursor) == 2 && col(*cursor) == 2 && value(*cursor) != 7) throw test_exception();
	if (row(*cursor) == 2 && col(*cursor) == 4 && value(*cursor) != 0) throw test_exception();
    }
}



template <typename Matrix, typename Tag, typename Complexity>
void two_d_iteration_impl(char const* outer, Matrix & matrix, Tag, Complexity)
{
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type                                 row(matrix); 
    typename traits::col<Matrix>::type                                 col(matrix); 
    typename traits::const_value<Matrix>::type                         value(matrix); 
    typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;
    // typedef typename traits::range_generator<Tag, Matrix>::complexity  complexity;

    cout << outer << ": " << Complexity() << '\n';
    // check_same_type(complexity(), ExpComplexity());
    for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	typedef mtl::tag::nz     inner_tag;
	cout << "---\n";
	typedef typename traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor)
	    cout << "matrix[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << '\n';
    }
} 


template <typename Matrix, typename Tag>
void two_d_iteration_impl(char const* name, Matrix &, Tag, mtl::complexity_classes::infinite)
{
    cout << name << ": Tag has no implementation\n";
}

template <typename Matrix, typename Tag, typename Complexity>
void two_d_iterator_iteration_impl(char const* outer, Matrix & matrix, Tag, Complexity)
{
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type                                 row(matrix); 
    typename traits::col<Matrix>::type                                 col(matrix); 
    typename traits::const_value<Matrix>::type                         value(matrix); 
    typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;

    cout << outer << " with iterators: " << Complexity() << '\n';
    for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	typedef mtl::tag::const_iter::nz     inner_tag;
	cout << "---\n";
	typedef typename traits::range_generator<inner_tag, cursor_type>::type iter_type;
	for (iter_type iter = mtl::begin<inner_tag>(cursor), i_end = mtl::end<inner_tag>(cursor); iter != i_end; ++iter)
	    cout << *iter << '\n';
    }
} 


template <typename Matrix, typename Tag>
void two_d_iterator_iteration_impl(char const* name, Matrix &, Tag, mtl::complexity_classes::infinite)
{
    cout << name << ": Tag has no implementation\n";
}

template <typename Matrix, typename Tag>
void two_d_iteration(char const* name, Matrix & matrix, Tag)
{
    typedef typename mtl::traits::range_generator<Tag, Matrix>::complexity  complexity;
    two_d_iteration_impl(name, matrix, Tag(), complexity());
    two_d_iterator_iteration_impl(name, matrix, Tag(), complexity());
}    




template <typename Matrix>
void matrix_init(Matrix& matrix)
{
    typedef typename Matrix::parameters   parameters;
    typedef typename Matrix::value_type   value_type;

    mtl::mat::compressed2D_inserter<value_type, parameters> inserter(matrix);
    inserter(2, 2) << 7; inserter(1, 4) << 3; inserter(3, 2) << 9; inserter(5, 1) << 5;
}
    
 
template <typename Orientation, typename Indexing>
void test_compressed2D(char const* name)
{
    cout << "\n====================\n" << name << "\n====================\n";
    typedef mtl::mat::parameters<Orientation, Indexing, mtl::fixed::dimensions<8, 6> >   parameters;
    typedef mtl::compressed2D<int, parameters>                                              matrix_type;
    matrix_type                                                                             matrix; 

    matrix_init(matrix);
    std::cout << "\n\n";
    print_matrix(matrix);

    one_d_iteration("\nMatrix", matrix); 
    //one_d_iterator_iteration("\nMatrix (iterator)", matrix); 

    two_d_iteration("Row-wise", matrix, mtl::tag::row());
    two_d_iteration("Column-wise", matrix, mtl::tag::col());
    two_d_iteration("On Major", matrix, mtl::tag::major());


    mtl::transposed_view<matrix_type> trans_matrix(matrix);
    cout << "\n===\n";
    print_matrix(trans_matrix);

    one_d_iteration("\nTransposed matrix", trans_matrix);
    //one_d_iterator_iteration("\nMatrix (iterator)", trans_matrix); 

    two_d_iteration("Transposed row-wise", trans_matrix, mtl::tag::row());
    two_d_iteration("Transposed Column-wise", trans_matrix, mtl::tag::col());
    two_d_iteration("Transposed On Major", trans_matrix, mtl::tag::major());
}

int main(int, char**)
{
    test_compressed2D<mtl::row_major, mtl::index::c_index>("CRS");
    // test_compressed2D<row_major, mtl::index::f_index>("CRS Fortran"); deprecated
    test_compressed2D<mtl::col_major, mtl::index::c_index>("CCS");
    // test_compressed2D<col_major, mtl::index::f_index>("CCS Fortran"); deprecated
    
    return 0;
}
 
