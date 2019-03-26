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
#include <string>
#include <boost/tuple/tuple.hpp>

#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp>
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/operation/print_matrix.hpp>
#include <boost/numeric/mtl/operation/sub_matrix.hpp>
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>
#include <boost/numeric/mtl/recursion/base_case_test.hpp>
#include <boost/numeric/mtl/recursion/for_each.hpp>



using namespace std;  


template <typename Recursator>
void print_depth_first(Recursator const& recursator, string str)
{
    if (recursator.is_empty())
	return;
    cout << "\nRecursion: " << str << endl << *recursator;

  
    // for full recursion remove the string length limitation
    if (!recursator.is_leaf()) { 
	print_depth_first(recursator.north_west(), string("north west of ") + str);
	print_depth_first(recursator.south_west(), string("south west of ") + str);
	print_depth_first(recursator.north_east(), string("north east of ") + str);
	print_depth_first(recursator.south_east(), string("south east of ") + str);
    }
} 


template <typename Recursator, typename BaseCaseTest>
void recursive_print(Recursator const& recursator, string str, BaseCaseTest const& is_base)
{
    if (is_base(recursator))
	cout << "\nBase case: " << str << endl << *recursator;
    else {
	recursive_print(recursator.north_west(), string("north west of ") + str, is_base);
	recursive_print(recursator.south_west(), string("south west of ") + str, is_base);
	recursive_print(recursator.north_east(), string("north east of ") + str, is_base);
	recursive_print(recursator.south_east(), string("south east of ") + str, is_base);
    }
} 


template <typename Recursator, typename BaseCaseTest>
void recursive_print_checked(Recursator const& recursator, string str, BaseCaseTest const& is_base)
{
    if (recursator.is_empty())
	return;
    if (is_base(recursator)) {
	cout << "\nBase case: " << str << endl << *recursator;
    } else {
	recursive_print_checked(recursator.north_west(), string("north west of ") + str, is_base);
	recursive_print_checked(recursator.south_west(), string("south west of ") + str, is_base);
	recursive_print_checked(recursator.north_east(), string("north east of ") + str, is_base);
	recursive_print_checked(recursator.south_east(), string("south east of ") + str, is_base);
    }
} 

struct print_functor
{
    template <typename Matrix>
    void operator() (Matrix const& matrix) const
    {
	cout << matrix << endl;
    }
};

template <typename Matrix>
void test_sub_matrix(Matrix& matrix)
{
    using mtl::recursion::for_each; using mtl::recursion::max_dim_test; 
    using mtl::mat::recursator; using mtl::transposed_view;

    cout << matrix << endl;    
    
    max_dim_test              is_base(2);
    recursator<Matrix>        rec(matrix);
    recursive_print_checked(rec, "", is_base);
	 
    cout << "\n====================\n"
	 <<   "Same with transposed\n"
	 <<   "====================\n\n";

    transposed_view<Matrix> trans_matrix(matrix);

    cout << trans_matrix; 
    recursator< transposed_view<Matrix> > trans_recursator(trans_matrix);
    // print_depth_first(trans_recursator, "");
    recursive_print_checked(trans_recursator, "", is_base);
	 
    cout << "\n=============================\n"
	 <<   "Again with recursive for_each\n"
	 <<   "=============================\n\n";

    for_each(trans_recursator, print_functor(), is_base);
}


template <typename Matrix>
void fill_matrix(Matrix& matrix)
{
    namespace traits = mtl::traits;
    using mtl::begin; using mtl::end;

    typename traits::row<Matrix>::type                                 row(matrix);
    typename traits::col<Matrix>::type                                 col(matrix);
    typename traits::value<Matrix>::type                               value(matrix);
    typedef  glas::tag::nz                                          tag;
    typedef typename traits::range_generator<tag, Matrix>::type        cursor_type;
    
    double x= 10.3;
    for (cursor_type cursor = begin<tag>(matrix), cend = end<tag>(matrix); cursor != cend; ++cursor) {
	value(*cursor, x);
	x+= 1.0; 
    }
       
}
  
 
int main(int, char**)
{
    using namespace mtl;

    cout << "=====================\n"
	 << "Morton-ordered matrix\n"
	 << "=====================\n\n";

    typedef morton_dense<double,  0x55555555> matrix_type;    
    matrix_type matrix(6, 5);   
    fill_matrix(matrix); 
    test_sub_matrix(matrix);

    cout << "\n=========================\n"
	 << "Doppler matrix (4x4 base)\n"
	 << "=========================\n\n";

    typedef morton_dense<double,  0x55555553> dmatrix_type;    
    dmatrix_type dmatrix(6, 5);   
    fill_matrix(dmatrix); 
    test_sub_matrix(dmatrix);

    cout << "\n======================\n"
	 << "Row-major dense matrix\n"
	 << "======================\n\n";

    dense2D<double, mat::parameters<> >   rmatrix(non_fixed::dimensions(6, 5));
    fill_matrix(rmatrix); 
    test_sub_matrix(rmatrix);
 
    cout << "=================================\n"
	 << "Vector-like morton-ordered matrix\n"
	 << "=================================\n\n";

    matrix_type vmatrix(17, 2);   
    fill_matrix(vmatrix); 
    test_sub_matrix(vmatrix);

    return 0;
}

