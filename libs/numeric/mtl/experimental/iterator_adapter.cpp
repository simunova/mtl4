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
#include <boost/test/minimal.hpp>
#include <boost/tuple/tuple.hpp>
#include <numeric>

#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/utility/range_generator.hpp>
#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/operation/raw_copy.hpp>
#include <boost/numeric/mtl/utilities/iterator_adapter.hpp>

using namespace mtl;
using namespace std;

using mtl::utilities::iterator_adapter;
using mtl::utilities::const_iterator_adapter;

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
	typename traits::row<Matrix>::type                         r = row(matrix);
	typename traits::col<Matrix>::type                         c = col(matrix);
	typename traits::value<Matrix>::type                       v = value(matrix);
	typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;
	typedef typename traits::range_generator<Tag, Matrix>::complexity  complexity;

	cout << outer << complexity() << '\n';
	check_same_type(complexity(), ExpComplexity());
	for (cursor_type cursor = begin<Tag>(matrix), cend = end<Tag>(matrix); cursor != cend; ++cursor) {
	    typedef glas::tag::all     inner_tag;
	    typedef typename traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	    for (icursor_type icursor = begin<inner_tag>(cursor), icend = end<inner_tag>(cursor); icursor != icend; ++icursor)
		cout << "matrix[" << r(*icursor) << ", " << c(*icursor) << "] = " << v(*icursor) << '\n';
	}
    }

    template <typename Matrix>
    void one_d_iteration(char const* name, Matrix & matrix, size_t check_row, size_t check_col, double check)
    {
	typename traits::row<Matrix>::type                         r = row(matrix);
	typename traits::col<Matrix>::type                         c = col(matrix);
	typename traits::value<Matrix>::type                       v = value(matrix);
	typedef  glas::tag::nz                                  tag;
	typedef typename traits::range_generator<tag, Matrix>::type        cursor_type;
	typedef typename traits::range_generator<tag, Matrix>::complexity  complexity;

	cout << name << "\nElements: " << complexity() << '\n';

	for (cursor_type cursor = begin<tag>(matrix), cend = end<tag>(matrix); cursor != cend; ++cursor) {
	    cout << "matrix[" << r(*cursor) << ", " << c(*cursor) << "] = " << v(*cursor) << '\n';
	    if (r(*cursor) == check_row && c(*cursor) == check_col && v(*cursor) != check) throw test_dense2D_exception();
	}

	typedef const_iterator_adapter<typename traits::row<Matrix>::type, cursor_type, double> row_adapter_type;
	row_adapter_type row_iter(r, begin<tag>(matrix)), row_end(r, end<tag>(matrix));
	typedef const_iterator_adapter<typename traits::col<Matrix>::type, cursor_type, double> col_adapter_type;
	col_adapter_type col_iter(c, begin<tag>(matrix));
	typedef iterator_adapter<typename traits::value<Matrix>::type, cursor_type, double> value_adapter_type;
	value_adapter_type value_iter(v, begin<tag>(matrix));

	cout << "\nSame with iterator adapter.\n";
	for (; row_iter != row_end; ++row_iter, ++col_iter, ++value_iter) {
	    cout << "matrix[" << *row_iter << ", " << *col_iter << "] = " << *value_iter << '\n';
	}

	value_adapter_type value_first(v, begin<tag>(matrix)), value_last(v, end<tag>(matrix));
	cout << "Sum of all values = " << accumulate(value_first, value_last, 0.0) << "\n";
    }
    
    void operator() (double element_1_2)
    {
	typedef dense2D<double, Parameters> matrix_type;
	matrix_type   matrix;
	double        val[] = {1., 2., 3., 4., 5., 6.};
	raw_copy(val, val+6, matrix);

	one_d_iteration("\nMatrix", matrix, 1, 2, element_1_2);
	two_d_iteration("\nRows: ", matrix, glas::tag::row(), ExpRowComplexity());
	two_d_iteration("\nColumns: ", matrix, glas::tag::col(), ExpColComplexity());

	transposed_view<matrix_type> trans_matrix(matrix);
	one_d_iteration("\nTransposed matrix", trans_matrix, 2, 1, element_1_2);
	two_d_iteration("\nRows: ", trans_matrix, glas::tag::row(), ExpColComplexity());
	two_d_iteration("\nColumns: ", trans_matrix, glas::tag::col(), ExpRowComplexity());
    }
};

int test_main(int argc, char* argv[])
{
    typedef mat::parameters<row_major, mtl::index::c_index, fixed::dimensions<2, 3> > parameters1;
    test_dense2D<parameters1, complexity::linear_cached, complexity::linear>()(6.0);

    typedef mat::parameters<row_major, mtl::index::f_index, fixed::dimensions<2, 3> > parameters2;
    test_dense2D<parameters2, complexity::linear_cached, complexity::linear>()(2.0);

    typedef mat::parameters<col_major, mtl::index::c_index, fixed::dimensions<2, 3> > parameters3;
    test_dense2D<parameters3, complexity::linear, complexity::linear_cached>()(6.0);

    typedef mat::parameters<col_major, mtl::index::f_index, fixed::dimensions<2, 3> > parameters4;
    test_dense2D<parameters4, complexity::linear, complexity::linear_cached>()(3.0);

    return 0;
}
