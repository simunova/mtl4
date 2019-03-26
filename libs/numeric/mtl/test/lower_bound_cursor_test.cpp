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

#include <boost/numeric/mtl/mtl.hpp>


using namespace std;

struct test_exception {};

// Might be specialized on Collection to replace < by <>

template <typename Cursor>
struct check_position_aux
{
    template <typename Matrix, typename Coll>
    void operator()(const Matrix& A, const Coll& c, const Cursor& cursor, bool check_row) const
    {
	namespace traits= mtl::traits;
	typename traits::row<Matrix>::type           row(A); 
	typename traits::col<Matrix>::type           col(A); 
	std::cout << ", cursor is pointing at [" << row(*cursor) << "][" << col(*cursor) << "]\n";

	typename traits::range_generator<mtl::tag::nz, Coll>::type bref= mtl::begin<mtl::tag::nz>(c), eref= mtl::end<mtl::tag::nz>(c);

	if (check_row) {
	    // same column and row at value (or beyond) or end cursor 
	    // if ( (!(col(*cursor) == col(*bref) && row(*cursor) >= 2) || cursor == eref ) )
	    MTL_THROW_IF(col(*cursor) != col(*bref) || !(cursor == eref || row(*cursor) >= 2), mtl::runtime_error("Cursor's row must be 2 (or larger)"));
	} else
	    // same row and column at value (or beyond) or end cursor
	    MTL_THROW_IF(row(*cursor) != row(*bref) || !(cursor == eref || col(*cursor) >= 2), mtl::runtime_error("Cursor's column must be 2 (or larger)"));
    }
};


template <typename Matrix, typename Cursor, int Complexity>
struct check_position_aux<mtl::traits::detail::sub_matrix_cursor<Matrix, Cursor, Complexity> >
{
    template <typename Matrix2, typename Coll>
    void operator()(const Matrix2&, const Coll&, const mtl::traits::detail::sub_matrix_cursor<Matrix, Cursor, Complexity>& cursor, bool ) const
    {
	std::cout << ", cursor is pointing at " << cursor.value() << "\n";
	MTL_THROW_IF(cursor.value() < 2, mtl::runtime_error("Cursor must be 2 (or larger)"));
    }
};


template <typename Matrix, typename Coll, typename Cursor>
void check_position(const Matrix& A, const Coll& c, const Cursor& cursor, bool check_row)
{
    check_position_aux<Cursor>()(A, c, cursor, check_row);
}


template <typename Matrix, typename Coll, typename Tag>
void test(char const* name, const Matrix& A, const Coll& c, Tag, bool cr)
{
    namespace traits= mtl::traits; using mtl::lower_bound;

    typename traits::range_generator<Tag, Coll>::type  cursor= lower_bound<Tag>(c, 2);

    //std::cout << "Type of cursor is " << typeid(cursor).name() << ", type of key is " << typeid(*cursor).name() << "\n";
    std::cout << name; 
    check_position(A, c, cursor, cr);
}



template <typename Matrix>
void rm_matrix_test(char const* name, const Matrix& matrix)
{
    using mtl::tag::row; using mtl::tag::nz; using mtl::traits::range_generator;
    std::cout << "\n\nTesting matrix type " << name << "\n";
    
    test("Find row in matrix", matrix, matrix, row(), true);

    typename range_generator<row, Matrix>::type row_cursor(mtl::begin<row>(matrix));
    test("Find non-zero in row 0", matrix, row_cursor, nz(), false);
    ++row_cursor; test("Find non-zero in row 1", matrix, row_cursor, nz(), false);
    ++row_cursor; test("Find non-zero in row 2", matrix, row_cursor, nz(), false);
}

template <typename Matrix>
void cm_matrix_test(char const* name, const Matrix& matrix)
{
    using mtl::tag::col; using mtl::tag::nz; using mtl::traits::range_generator;
    std::cout << "\n\nTesting matrix type " << name << "\n";
    
    test("Find column in matrix", matrix, matrix, col(), false);

    typename range_generator<col, Matrix>::type col_cursor(mtl::begin<col>(matrix));
    test("Find non-zero in column 0", matrix, col_cursor, nz(), true);
    ++col_cursor; test("Find non-zero in column 1", matrix, col_cursor, nz(), true);
    ++col_cursor; test("Find non-zero in column 2", matrix, col_cursor, nz(), true);
}

template <typename Matrix>
void dense_matrix_test(char const* name, const Matrix& matrix)
{
    using mtl::tag::row; using mtl::tag::col; using mtl::tag::nz; using mtl::traits::range_generator;
    std::cout << "\n\nTesting matrix type " << name << "\n";
    
    test("Find row in matrix", matrix, matrix, row(), true);

    typename range_generator<row, Matrix>::type row_cursor(mtl::begin<row>(matrix));
    test("Find non-zero in row 0", matrix, row_cursor, nz(), false);
    ++row_cursor; test("Find non-zero in row 1", matrix, row_cursor, nz(), false);
    ++row_cursor; test("Find non-zero in row 2", matrix, row_cursor, nz(), false);

    test("Find column in matrix", matrix, matrix, col(), false);

    typename range_generator<col, Matrix>::type col_cursor(mtl::begin<col>(matrix));
    test("Find non-zero in column 0", matrix, col_cursor, nz(), true);
    ++col_cursor; test("Find non-zero in column 1", matrix, col_cursor, nz(), true);
    ++col_cursor; test("Find non-zero in column 2", matrix, col_cursor, nz(), true);
}



int main(int, char**)
{
    using namespace mtl;
    typedef mat::parameters<col_major> col_para;

    double ar[4][4]= {{1., 0., 2., 0.},
		     {0., 3., 0., 4.},
		     {0., 5., 0., 0.},
		     {0., 0., 6., 0.}};

    compressed2D<double>                cr(ar); crop(cr);
    compressed2D<double, col_para>      cc(ar); crop(cc);       
    dense2D<double>                     dr(ar);
    dense2D<double, col_para>           dc(ar);
    morton_dense<double, recursion::morton_z_mask> mo(ar);

    rm_matrix_test("compressed2D<double>",                   cr);
    cm_matrix_test("compressed2D<double, col_para>",         cc); 
    dense_matrix_test("dense2D<double>",                     dr);
    dense_matrix_test("dense2D<double, col_para>",           dc);
    dense_matrix_test("morton_dense<double, recursion::morton_z_mask>", mo);

    return 0;
}
 
