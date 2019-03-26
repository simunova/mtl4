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
#include <boost/tuple/tuple.hpp>

#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/operation/raw_copy.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>



using namespace std;  


struct test_morton_dense 
{
    template <typename Matrix, typename Tag>
    void two_d_iteration(char const* outer, Matrix & matrix, Tag)
    {
	typename mtl::traits::row<Matrix>::type                                 row(matrix);
	typename mtl::traits::col<Matrix>::type                                 col(matrix);
	typename mtl::traits::value<Matrix>::type                               value(matrix);
	typedef typename mtl::traits::range_generator<Tag, Matrix>::type        cursor_type;

	cout << outer << '\n';
	for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	    typedef mtl::tag::all     inner_tag;
	    typedef typename mtl::traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	    for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor)
		cout << "matrix[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << '\n';
	    icursor_type ibeg = mtl::begin<inner_tag>(cursor), icursor= ibeg + 2;
	    cout << "--\nmatrix[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << "\n--\n";
	}
    }

    template <typename Matrix, typename Tag>
    void two_d_iterator_iteration(char const* outer, Matrix & matrix, Tag)
    {
	typename mtl::traits::row<Matrix>::type                                 row(matrix);
	typename mtl::traits::col<Matrix>::type                                 col(matrix);
	typename mtl::traits::value<Matrix>::type                               value(matrix);
	typedef typename mtl::traits::range_generator<Tag, Matrix>::type        cursor_type;

	cout << outer << '\n';
	for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	    typedef mtl::tag::iter::all     inner_tag;
	    typedef typename mtl::traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	    for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor)
		cout << *icursor <<'\n';
	}
    } 

    template <typename Matrix> 
    void one_d_iteration(char const* name, Matrix & matrix)
    {
	typename mtl::traits::row<Matrix>::type                                 row(matrix);
	typename mtl::traits::col<Matrix>::type                                 col(matrix);
	typename mtl::traits::value<Matrix>::type                               value(matrix);
	typedef  mtl::tag::nz                                          tag; 
	typedef typename mtl::traits::range_generator<tag, Matrix>::type        cursor_type;

	cout << name << "\nElements: \n";
	for (cursor_type cursor = mtl::begin<tag>(matrix), cend = mtl::end<tag>(matrix); cursor != cend; ++cursor) {
	    cout << "matrix[" << row(*cursor) << ", " << col(*cursor) << "] = " << value(*cursor) << '\n';
	}
    }
    
    template <typename Matrix>
    void fill_matrix(Matrix & matrix)
    {
	typename mtl::traits::value<Matrix>::type                               value(matrix);
	typedef  mtl::tag::nz                                          tag;
	typedef typename mtl::traits::range_generator<tag, Matrix>::type        cursor_type;

	typename Matrix::value_type  v= 1;

	for (cursor_type cursor = mtl::begin<tag>(matrix), cend = mtl::end<tag>(matrix); cursor != cend; ++cursor) {
	    value(*cursor, v);
	    v+= 1;
	}
    }

    template <typename Matrix>
    void check_cursor_increment(Matrix& matrix)
    {
	typename mtl::traits::row<Matrix>::type                                 row(matrix);
	typename mtl::traits::col<Matrix>::type                                 col(matrix);
	typename mtl::traits::value<Matrix>::type                               value(matrix);
	typedef  mtl::tag::nz                                          tag;
	typedef typename mtl::traits::range_generator<tag, Matrix>::type        cursor_type;
	
	cursor_type cursor = mtl::begin<tag>(matrix);
	cout << "begin: matrix[" << row(*cursor) << ", " << col(*cursor) << "] = " << value(*cursor) << '\n';
	cursor.advance(2, 2);
	cout << "advance (2,2): matrix[" << row(*cursor) << ", " << col(*cursor) << "] = " << value(*cursor) << '\n';
	cursor.advance(-1, -1);
	cout << "advance (-1, -1): matrix[" << row(*cursor) << ", " << col(*cursor) << "] = " << value(*cursor) << '\n';
    }

    template <typename Matrix>
    void operator() (Matrix& matrix)
    {
	fill_matrix(matrix);
	check_cursor_increment(matrix);

	one_d_iteration("\nMatrix", matrix);
	two_d_iteration("\nRows: ", matrix, mtl::tag::row());
	two_d_iteration("\nColumns: ", matrix, mtl::tag::col());
	two_d_iterator_iteration("\nRows (iterator): ", matrix, mtl::tag::row());
	two_d_iterator_iteration("\nColumns (iterator): ", matrix, mtl::tag::col());

	mtl::transposed_view<Matrix> trans_matrix(matrix);
	one_d_iteration("\nTransposed matrix", trans_matrix);
	two_d_iteration("\nRows: ", trans_matrix, mtl::tag::row());
	two_d_iteration("\nColumns: ", trans_matrix, mtl::tag::col());
	two_d_iterator_iteration("\nRows (iterator): ", trans_matrix, mtl::tag::row());
	two_d_iterator_iteration("\nColumns (iterator): ", trans_matrix, mtl::tag::col());
    }
};



 
int main(int, char**)
{
    using namespace mtl;

    morton_dense<double,  0x55555555> matrix1(3, 5);
    matrix1[1][3]= 2.3;
    cout << "matrix1[1][3] = " << matrix1[1][3] << endl;

    typedef morton_dense<int,  0x55555553> matrix2_type;
    matrix2_type                           matrix2(5, 6);
    matrix2[1][3]= 3;
    cout << "matrix2[1][3] = " << matrix2[1][3] << endl;
    

    test_morton_dense()(matrix1);
    test_morton_dense()(matrix2);

    return 0;
}
