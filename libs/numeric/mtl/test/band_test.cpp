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

// #define MTL_DEBUG_DMAT_DMAT_MULT // For debugging evidently

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;

template <typename Matrix>
void list_entries(const Matrix& A, int begin, int end)
{
    using namespace mtl; using namespace mtl::tag; using mtl::tag::major; using mtl::traits::range_generator;  
    typedef typename range_generator<major, Matrix>::type     cur_type;    
    typedef typename range_generator<nz, cur_type>::type      icur_type; 
    typename mtl::traits::col<Matrix>::type                   col(A);
    typename mtl::traits::row<Matrix>::type                   row(A);
    typename mtl::traits::const_value<Matrix>::type           value(A);

    for(cur_type c= mtl::begin<major>(A), cend= mtl::end<major>(A); c != cend; ++c)
	for(icur_type kc= mtl::begin<nz>(c), kend= mtl::end<nz>(c); kc != kend; ++kc) {
	    std::cout << "A[" << row(*kc) << "][" << col(*kc) << "] = " << value(*kc) << "\n";
	    int band= col(*kc) - row(*kc);
	    if (band < begin) std::cout << "outside on the left!\n";
	    if (band >= end) std::cout << "outside on the right!\n";
	}
}

template <typename Matrix>
void check(const Matrix& A, int begin, int end)
{
    // list_entries(A, begin, end);
    typedef typename mtl::Collection<Matrix>::value_type   value_type;
    typedef typename mtl::Collection<Matrix>::size_type    size_type;

    for (size_type i= 0; i < num_rows(A); i++)
	for (size_type j= 0; j < num_cols(A); j++) {
	    long band= long(j) - long(i);
	    if (band < begin) {
		MTL_THROW_IF(A[i][j] != value_type(0), mtl::runtime_error("Value must be zero left of the bands"));
	    } else if (band >= end) {
		MTL_THROW_IF(A[i][j] != value_type(0), mtl::runtime_error("Value must be zero right of the bands"));
	    } else
		MTL_THROW_IF(A[i][j] != value_type(band + 0), mtl::runtime_error("Wrong non-zero value within the bands"));
	}
}


template <typename Matrix>
void test(Matrix& A, const char* name)
{
    using mtl::mat::inserter; using mtl::Collection;

    typedef typename Collection<Matrix>::value_type   value_type;

    A.change_dim(6, 5);
    {
	inserter<Matrix>   ins(A);
	for (unsigned i= 0; i < num_rows(A); i++)
	    for (unsigned j= 0; j < num_cols(A); j++)
		ins[i][j]= value_type(j) - value_type(i);
    }
    cout << "\n" << name << "\n";  
    cout << "\n" << name << "\n" << "A =\n" << A;
    Matrix B( bands(A, 2, 4) );
    cout << "\nbands(A, 2, 4) = \n" << B;
    check(B, 2, 4);

    Matrix U( upper(A) );
    cout << "\nupper(A) = \n" << U;
    check(U, 0, 10000);

    Matrix SU( strict_upper(A) );
    cout << "\nstrict_upper(A) = \n" << SU;
    check(SU, 1, 10000);

    Matrix L( lower(A) );
    cout << "\nlower(A) = \n" << L;
    check(L, -10000, 1);
    
    Matrix SL( strict_lower(A) );
    cout << "\nstrict_lower(A) = \n" << SL;
    check(SL, -10000, 0);
 
    Matrix P( trans(A) * upper(A) );
    cout << " trans(A) * upper(A) = \n" << with_format(P, 4, 3);
    Matrix P_cmp( trans(A) * U );
    Matrix P_cmp2( trans(A) * A );
    cout << "\ntrans(A) * upper(A) = \n" << with_format(P, 4, 3);
    cout << " for comparison trans(A) * U = \n" << with_format(P_cmp, 4, 3);
    cout << " for comparison trans(A) * A = \n" << with_format(P_cmp2, 4, 3);
    MTL_THROW_IF(abs(P[1][1] - P_cmp[1][1]) > .00001, mtl::runtime_error("Multiplication wrong"));

#if 0
    // Take this out as sparse matrices have no sub-matrix; only for debugging wrong products
    Matrix A2= sub_matrix(A, 0, 3, 0, 3), U2= upper(A2);
    Matrix P2= A2 * upper(A2); //, P2_cmp= A2 * U2, A2_square= A2 * A2;
    cout << "\nA2 * upper(A2) = \n" << with_format(P2, 4, 3);
    cout << " for comparison A2 * U = \n" << with_format(P2_cmp, 4, 3);
    cout << " for comparison A2 * A2 = \n" << with_format(A2_square, 4, 3);
    MTL_THROW_IF(abs(P2[1][1] - P2_cmp[1][1]) > .00001, mtl::runtime_error("Multiplication wrong"));
#endif


#if 0
    // Would too painfully slow !

    Matrix D( diagonal(A) );
    cout << "\ndiagonal(A) = \n" << D;
    check(D, 0, 1);

    Matrix T( tri_diagonal(A) );
    cout << "\ntri_diagonal(A) = \n" << T;
    check(T, 0, 1);
#endif
}


int main(int, char**)
{
    using namespace mtl;

    dense2D<double>                                      dr;
    dense2D<double, mat::parameters<col_major> >      dc;
    morton_dense<double, recursion::morton_z_mask>       mzd;
    morton_dense<double, recursion::doppled_2_row_mask>  d2r;
    compressed2D<double>                                 cr;
    compressed2D<double, mat::parameters<col_major> > cc;

    dense2D<complex<double> >                            drc;
    compressed2D<complex<double> >                       crc;

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(drc, "Dense row major complex");
    test(crc, "Compressed row major complex");

    // For better readability I don't want finish with a complex
    test(cc, "Compressed column major");

    return 0;
}
