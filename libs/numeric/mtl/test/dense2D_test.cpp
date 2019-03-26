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

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/operation/raw_copy.hpp>
#include <boost/numeric/mtl/operation/sub_matrix.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/matrix/inserter.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>

// This test is old and ugly, consider replacing it without decreasing coverage

using namespace std;

struct test_dense2D_exception {};

template <typename T1, typename T2>
void check_same_type(T1, T2)
{
    throw test_dense2D_exception();
}
 
// If same type we're fine
template <typename T1>
void check_same_type(T1, T1) {}


template <typename Parameters, typename ExpRowComplexity, typename ExpColComplexity>
struct test_dense2D
{
    template <typename Matrix, typename Tag, typename ExpComplexity>
    void two_d_iteration(char const* outer, Matrix & matrix, Tag, ExpComplexity)
    {
	namespace traits = mtl::traits;
	typename traits::row<Matrix>::type                                 row(matrix); 
	typename traits::col<Matrix>::type                                 col(matrix); 
	typename traits::const_value<Matrix>::type                         value(matrix); 
	typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;
	typedef typename traits::range_generator<Tag, Matrix>::complexity  complexity;

	cout << outer << complexity() << '\n';
	check_same_type(complexity(), ExpComplexity());
	for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	    typedef mtl::tag::all     inner_tag;
	    typedef typename traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	    for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor)
		cout << "matrix[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << '\n';
	}
    } 

    template <typename Matrix, typename Tag, typename ExpComplexity>
    void two_d_iterator_iteration(char const* outer, Matrix & matrix, Tag, ExpComplexity)
    {
	namespace traits = mtl::traits;
	typename traits::row<Matrix>::type                                 row(matrix); 
	typename traits::col<Matrix>::type                                 col(matrix); 
	typename traits::const_value<Matrix>::type                         value(matrix); 
	typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;
	typedef typename traits::range_generator<Tag, Matrix>::complexity  complexity;

	cout // << "Matrix traversal with iterators" 
	     << outer << complexity() << '\n';
	check_same_type(complexity(), ExpComplexity());
	for (cursor_type cursor = mtl::begin<Tag>(matrix), cend = mtl::end<Tag>(matrix); cursor != cend; ++cursor) {
	    typedef mtl::tag::iter::all     inner_tag;
	    typedef typename traits::range_generator<inner_tag, cursor_type>::type iter_type;
	    for (iter_type iter = mtl::begin<inner_tag>(cursor), i_end = mtl::end<inner_tag>(cursor); iter != i_end; ++iter)
		cout << *iter << '\n';
	}
    }


    template <typename Matrix>
    void one_d_iteration(char const* name, Matrix & matrix, size_t check_row, size_t check_col, double check)
    {
	namespace traits = mtl::traits;
	typename traits::row<Matrix>::type                                 row(matrix);
	typename traits::col<Matrix>::type                                 col(matrix);
	typename traits::value<Matrix>::type                               value(matrix); 
	typedef  mtl::tag::nz                                          tag;
	typedef typename traits::range_generator<tag, Matrix>::type        cursor_type;
	typedef typename traits::range_generator<tag, Matrix>::complexity  complexity;

	cout << name << "\nElements: " << complexity() << '\n';
	for (cursor_type cursor = mtl::begin<tag>(matrix), cend = mtl::end<tag>(matrix); cursor != cend; ++cursor) {
	    cout << "matrix[" << row(*cursor) << ", " << col(*cursor) << "] = " << value(*cursor) << '\n';
	    if (row(*cursor) == check_row && col(*cursor) == check_col && value(*cursor) != check) throw test_dense2D_exception();
	}
    }
    
    void operator() (double element_1_2)
    {
	typedef mtl::dense2D<double, Parameters> matrix_type;
	matrix_type   matrix;
	matrix= 1., 2., 3., 
	        4., 5., 6.;
	// double        val[] = {1., 2., 3., 4., 5., 6.};
	// mtl::raw_copy(val, val+6, matrix);
 
	one_d_iteration("\nMatrix", matrix, 1, 2, element_1_2);
	two_d_iteration("\nRows: ", matrix, mtl::tag::row(), ExpRowComplexity());
	two_d_iteration("\nColumns: ", matrix, mtl::tag::col(), ExpColComplexity());
	two_d_iterator_iteration("\nRows (iterator): ", matrix, mtl::tag::row(), ExpRowComplexity());
	two_d_iterator_iteration("\nColumns (iterator): ", matrix, mtl::tag::col(), ExpColComplexity());

	mtl::transposed_view<matrix_type> trans_matrix(matrix); 
	one_d_iteration("\nTransposed matrix", trans_matrix, 2, 1, element_1_2);
	two_d_iteration("\nRows: ", trans_matrix, mtl::tag::row(), ExpColComplexity());
	two_d_iteration("\nColumns: ", trans_matrix, mtl::tag::col(), ExpRowComplexity());
	two_d_iterator_iteration("\nRows (iterator): ", trans_matrix, mtl::tag::row(), ExpColComplexity());
	two_d_iterator_iteration("\nColumns (iterator): ", trans_matrix, mtl::tag::col(), ExpRowComplexity());

	cout << "\nmatrix[1][2] = " << matrix[1][2] << "\n";
	matrix[1][2]= 18.0;
	cout << "matrix[1][2] = " << matrix[1][2] << "\n";
	cout << "trans_matrix[2][1] = " << trans_matrix[2][1] << "\n"; 
       

	mtl::mat::inserter<matrix_type>  i(matrix);
	i(1, 2) << 17.0;
	cout << "matrix[1, 2] = " << matrix(1, 2) << "\n";	 
    }
}; 

int main(int, char**)
{
    using namespace mtl;

    typedef mat::parameters<row_major, mtl::index::c_index, fixed::dimensions<2, 3> > parameters1;
    test_dense2D<parameters1, complexity_classes::linear_cached, complexity_classes::linear>()(6.0);

    // Don't use this Fortran indexing, it is absolutely incomplete
    // typedef mat::parameters<row_major, mtl::index::f_index, fixed::dimensions<2, 3> > parameters2;
    // test_dense2D<parameters2, complexity_classes::linear_cached, complexity_classes::linear>()(2.0);

    typedef mat::parameters<col_major, mtl::index::c_index, fixed::dimensions<2, 3> > parameters3;
    test_dense2D<parameters3, complexity_classes::linear, complexity_classes::linear_cached>()(6.0);

    // typedef mat::parameters<col_major, mtl::index::f_index, fixed::dimensions<2, 3> > parameters4;
    // test_dense2D<parameters4, complexity_classes::linear, complexity_classes::linear_cached>()(3.0);

#ifdef MTL_WITH_INITLIST
    dense2D<int> A = {{1, 2}, {3, 4}};
    std::stringstream ss;
    ss << A;
    
    std::cout << ss.str() << std::endl;
    std::string refsting = "[1 2]\n[3 4]\n";
    
    std::cout << "Equal to refstring: " << std::boolalpha << (refsting == ss.str()) << std::endl;
#endif    
    
    return 0;
}
