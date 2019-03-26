// Software License for MTL
// 
// Copyright (C) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#include <iostream>
#include <cmath>

// #define MTL_HAS_BLAS
// #define MTL_USE_OPTERON_OPTIMIZATION

#include <boost/numeric/mtl/utility/glas_tag.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/morton_dense.hpp> 
#include <boost/numeric/mtl/matrix/transposed_view.hpp>
#include <boost/numeric/mtl/recursion/bit_masking.hpp>
#include <boost/numeric/mtl/operation/print.hpp>
#include <boost/numeric/mtl/operation/dmat_dmat_mult.hpp>
#include <boost/numeric/mtl/operation/mult.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
#include <boost/numeric/mtl/matrix/hessian_setup.hpp>
#include <boost/numeric/mtl/operation/assign_mode.hpp>
#include <boost/numeric/mtl/operation/mult_assign_mode.hpp>
#include <boost/numeric/mtl/recursion/base_case_test.hpp>


using namespace std;  


template <typename MatrixA, typename MatrixB, typename MatrixC>
void test(MatrixA& A, MatrixB& B, MatrixC& C, const char* name)
{
    using mtl::assign::assign_sum; using mtl::assign::plus_sum; 
    using mtl::assign::minus_sum; using mtl::assign::mult_assign_mode; 
    using mtl::recursion::bound_test_static; using mtl::gen_recursive_dmat_dmat_mult_t;

    hessian_setup(A, 1.0);
    hessian_setup(B, 2.0);

    std::cout << "\n" << name << "  --- calling simple mult:\n"; std::cout.flush();
    typedef mtl::gen_dmat_dmat_mult_t<>  mult_t;
    mult_t                               mult;

    mult(A, B, C);
    // cout << "correct result is:\n" << with_format(C, 5, 3);
    check_hessian_matrix_product(C, A.num_cols());

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    typedef mtl::gen_dmat_dmat_mult_t<plus_sum>  add_mult_t;
    add_mult_t add_mult;

    add_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    typedef mtl::gen_dmat_dmat_mult_t<minus_sum>  minus_mult_t;
    minus_mult_t minus_mult;

    minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

#if 0
    std::cout << "\n" << name << "  --- calling mult with cursors and property maps:\n"; std::cout.flush();
    gen_cursor_dmat_dmat_mult_t<>  cursor_mult;

    cursor_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols());

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    gen_cursor_dmat_dmat_mult_t<plus_sum>  cursor_add_mult;

    cursor_add_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    gen_cursor_dmat_dmat_mult_t<minus_sum>  cursor_minus_mult; 

    cursor_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);
#endif 
    std::cout << "\n" << name << "  --- calling mult with tiling:\n"; std::cout.flush();
    typedef mtl::gen_tiling_dmat_dmat_mult_t<2, 2>  tiling_mult_t;
    tiling_mult_t tiling_mult;

    tiling_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols()); 

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    typedef mtl::gen_tiling_dmat_dmat_mult_t<2, 2, plus_sum>  tiling_add_mult_t;
    tiling_add_mult_t tiling_add_mult;

    tiling_add_mult(A, B, C); 
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    typedef mtl::gen_tiling_dmat_dmat_mult_t<2, 2, minus_sum>  tiling_minus_mult_t;
    tiling_minus_mult_t tiling_minus_mult;

    tiling_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

 
    std::cout << "\n" << name << "  --- calling mult with tiling 2x2:\n"; std::cout.flush();
    typedef mtl::gen_tiling_22_dmat_dmat_mult_t<>  tiling_22_mult_t;
    tiling_22_mult_t tiling_22_mult;

    tiling_22_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols()); 

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    // typedef gen_tiling_22_dmat_dmat_mult_t<plus_sum>  tiling_22_add_mult_t;
    typedef typename mult_assign_mode<tiling_22_mult_t, plus_sum>::type   tiling_22_add_mult_t;
    tiling_22_add_mult_t tiling_22_add_mult;

    tiling_22_add_mult(A, B, C); 
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    typedef mtl::gen_tiling_22_dmat_dmat_mult_t<minus_sum>  tiling_22_minus_mult_t;
    tiling_22_minus_mult_t tiling_22_minus_mult;

    tiling_22_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);


    std::cout << "\n" << name << "  --- calling mult with tiling 4x4:\n"; std::cout.flush();
    typedef mtl::gen_tiling_44_dmat_dmat_mult_t<>  tiling_44_mult_t;
    tiling_44_mult_t tiling_44_mult;

    tiling_44_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols()); 

    MatrixA A9(9, 9); MatrixB B9(9, 9); MatrixC C9(9, 9); // to test all block in 4x4 blocking for better coverage
    hessian_setup(A9, 1.0); hessian_setup(B9, 2.0);

    tiling_44_mult(A9, B9, C9);
    check_hessian_matrix_product(C9, num_cols(A9)); 

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    typedef mtl::gen_tiling_44_dmat_dmat_mult_t<plus_sum>  tiling_44_add_mult_t;
    tiling_44_add_mult_t tiling_44_add_mult;

    tiling_44_add_mult(A, B, C); 
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    typedef mtl::gen_tiling_44_dmat_dmat_mult_t<minus_sum>  tiling_44_minus_mult_t;
    tiling_44_minus_mult_t tiling_44_minus_mult;

    tiling_44_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

    typedef mtl::gen_recursive_dmat_dmat_mult_t<add_mult_t, bound_test_static<2>, plus_sum>  recursive_add_mult_t;

    std::cout << "\n" << name << "  --- calling mult recursively:\n"; std::cout.flush();
    // The recursive functor is C= A*B but the base case must be C+= A*B !!!!!!
    // gen_recursive_dmat_dmat_mult_t<add_mult_t, bound_test_static<32> >  recursive_mult;
    typename mult_assign_mode<recursive_add_mult_t, assign_sum>::type	recursive_mult;

    recursive_mult(A, B, C);
    // cout << "recursive result is:\n" << with_format(C, 5, 3);
    check_hessian_matrix_product(C, A.num_cols()); 

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    recursive_add_mult_t   recursive_add_mult;

    recursive_add_mult(A, B, C); 
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    // gen_recursive_dmat_dmat_mult_t<minus_mult_t, bound_test_static<32>, minus_sum>  recursive_minus_mult; 
    // check assign mode substitution both on matrix and on block level
    typename mult_assign_mode<recursive_add_mult_t, minus_sum>::type	recursive_minus_mult;
    recursive_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

    std::cout << "\n" << name << "  --- calling mult recursively with tiling:\n"; std::cout.flush();
    // The recursive functor is C= A*B but the base case must be C+= A*B !!!!!!
    gen_recursive_dmat_dmat_mult_t<tiling_add_mult_t, bound_test_static<32> >  recursive_tiling_mult;

    recursive_tiling_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols()); 
 
    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    gen_recursive_dmat_dmat_mult_t<tiling_add_mult_t, bound_test_static<32>, plus_sum>  recursive_tiling_add_mult;

    recursive_tiling_add_mult(A, B, C); 
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    gen_recursive_dmat_dmat_mult_t<tiling_minus_mult_t, bound_test_static<32>, minus_sum>  recursive_tiling_minus_mult; 

    recursive_tiling_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

    std::cout << "\n" << name << "  --- calling mult recursively platform specific plus tiling:\n"; std::cout.flush();
    typedef mtl::gen_platform_dmat_dmat_mult_t<plus_sum, tiling_add_mult_t> platform_tiling_add_mult_t;
    gen_recursive_dmat_dmat_mult_t<platform_tiling_add_mult_t, bound_test_static<32> >  recursive_platform_tiling_mult;
    
    recursive_platform_tiling_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

#ifdef MTL_HAS_BLAS
    std::cout << "\n" << name << "  --- calling blas mult:\n"; std::cout.flush(); 
    gen_blas_dmat_dmat_mult_t<>  blas_mult;
    blas_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols()); 
#endif    

#ifdef MTL_USE_OPTERON_OPTIMIZATION
    std::cout << "\n" << name << "  --- calling platform specific mult (empty):\n"; std::cout.flush(); 
    gen_platform_dmat_dmat_mult_t<>  platform_mult;
    platform_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols());

    std::cout << "\n" << name << "  --- check += :\n"; std::cout.flush();
    gen_platform_dmat_dmat_mult_t<plus_sum>  platform_add_mult;

    platform_add_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 2.0);
    
    std::cout << "\n" << name << "  --- check -= :\n"; std::cout.flush();
    gen_platform_dmat_dmat_mult_t<minus_sum>  platform_minus_mult;

    platform_minus_mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

#endif

    std::cout << "\n" << name << "  --- using mult(A, B, C) :\n"; std::cout.flush();
    mult(A, B, C);
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

    std::cout << "\n" << name << "  --- called as C= A * B:\n"; std::cout.flush();
    C= A * B;
    check_hessian_matrix_product(C, A.num_cols());

    std::cout << "\n" << name << "  --- check C+= A * B:\n"; std::cout.flush();
    C+= A * B;
    check_hessian_matrix_product(C, A.num_cols(), 2.0);

    std::cout << "\n" << name << "  --- check C-= A * B:\n"; std::cout.flush();
    C-= A * B;
    check_hessian_matrix_product(C, A.num_cols(), 1.0);

    if (A.num_cols() <= 10) 
	std::cout << A << "\n" << B << "\n" << with_format(C, 4, 4) << "\n";

}
 



template <typename MatrixA, typename MatrixB, typename MatrixC>
void single_test(MatrixA& A, MatrixB& B, MatrixC& C, const char*)
{
    using mtl::assign::plus_sum; using mtl::assign::minus_sum; 
    using mtl::recursion::bound_test_static;

    std::cout << "\n\n before matrix multiplication:\n";
    std::cout << "A:\n" << A;
    std::cout << "B:\n" << B;
    std::cout << "C:\n" << C << '\n';

    typedef mtl::gen_tiling_dmat_dmat_mult_t<2, 2, plus_sum>  tiling_add_mult_t;
    tiling_add_mult_t tiling_add_mult;
    tiling_add_mult(A, B, C); 
    
    std::cout << "\n\n after matrix multiplication:\n";
    std::cout << "A:\n" << A;
    std::cout << "B:\n" << B;
    std::cout << "C:\n" << C << '\n';
}

#ifdef MTL_HAS_BLAS
extern "C" {
void dgemm_(const char* transa, const char* transb, 
	    const int* m, const int* n, const int* k,
	    const double* alpha,  const double *da,  const int* lda,
	    const double *db, const int* ldb, const double* dbeta,
	    double *dc, const int* ldc);
}

typedef mtl::dense2D<double, mtl::mat::parameters<col_major> >        dc_t;

struct dgemm_t
{
    void operator()(const dc_t& A, const dc_t& B, dc_t& C)
    {
	int size= A.num_rows();
	double alpha= 1.0, beta= 0.0;
	dgemm_("N", "N", &size, &size, &size, &alpha, 
	       const_cast<double*>(&A[0][0]), &size, const_cast<double*>(&B[0][0]), 
	       &size, &beta, &C[0][0], &size);

    }
};


void test_blas()
{
    mtl::dense2D<double, mtl::mat::parameters<col_major> > A(7, 7), B(7, 7), C(7, 7);
    hessian_setup(A, 1.0);
    hessian_setup(B, 2.0);
    dgemm_t()(A, B, C);

    std::cout << C; 
    check_hessian_matrix_product(C, A.num_cols());
    
}

#endif // MTL_HAS_BLAS

int main(int argc, char* argv[])
{
    using namespace mtl;

    // Bitmasks:
    const unsigned long morton_mask= generate_mask<true, 0, row_major, 0>::value,
	doppled_32_row_mask_no_shark= generate_mask<true, 5, row_major, 0>::value,
	doppled_32_col_mask_no_shark= generate_mask<true, 5, col_major, 0>::value,
	doppled_32_row_mask= generate_mask<true, 5, row_major, 1>::value,
	doppled_32_col_mask= generate_mask<true, 5, col_major, 1>::value,
	doppled_z_32_row_mask= generate_mask<false, 5, row_major, 1>::value,
	doppled_z_32_col_mask= generate_mask<false, 5, col_major, 1>::value;
 
    unsigned size= 5; 
    if (argc > 1) size= atoi(argv[1]); 
    if (size < 2) size= 2;

    dense2D<double>                                  da(size, size-1), db(size-1, size-2), dc(size, size-2); 
    dense2D<double, mat::parameters<col_major> >  dca(size, size-1), dcb(size-1, size-2), dcc(size, size-2);
    dense2D<float>                                   fa(size, size-1), fb(size-1, size-2), fc(size, size-2);
    dense2D<float, mat::parameters<col_major> >   fca(size, size-1), fcb(size-1, size-2), fcc(size, size-2);
    morton_dense<double,  morton_mask>               mda(size, size-1), mdb(size-1, size-2), mdc(size, size-2);

    typedef morton_dense<double, doppled_32_row_mask_no_shark>  morton_t;
    morton_dense<double, doppled_32_row_mask_no_shark>      mrans(size, size-1), mrbns(size-1, size-2), mrcns(size, size-2);;
    morton_dense<double, doppled_32_col_mask_no_shark>      mcans(size, size-1), mcbns(size-1, size-2), mccns(size, size-2); 
    morton_dense<double, doppled_32_col_mask>      mca(size, size-1), mcb(size-1, size-2), mcc(size, size-2);
    morton_dense<double, doppled_32_row_mask>      mra(size, size-1), mrb(size-1, size-2), mrc(size, size-2);
    morton_dense<double, doppled_z_32_col_mask>    mzca(size, size-1), mzcb(size-1, size-2), mzcc(size, size-2);
    morton_dense<double, doppled_z_32_row_mask>    mzra(size, size-1), mzrb(size-1, size-2), mzrc(size, size-2);
    morton_dense<float, doppled_32_col_mask>       mcaf(size, size-1), mcbf(size-1, size-2), mccf(size, size-2);
    morton_dense<float, doppled_32_row_mask>       mraf(size, size-1), mrbf(size-1, size-2), mrcf(size, size-2);

    transposed_view<dense2D<double> > trans_db(db); 
    transposed_view<morton_t >        trans_mrbns(mrbns); 


    std::cout << "Testing different products\n";

#if 0
    test(da, trans_db, dc, "dense2D and transposed dense2D");
    test(mrans, trans_mrbns, mrcns, "hybrid with transposed matrix");
#endif

    test(da, db, dc, "dense2D");
    test(dca, dcb, dcc, "dense2D col-major");
    test(da, dcb, dc, "dense2D row x column-major");
    test(fa, fcb, fc, "dense2D float, row x column-major");
    test(da, fcb, fc, "dense2D mixed, dense and float"); 
    test(mda, mdb, mdc, "pure Morton");
    test(mca, mcb, mcc, "Hybrid col-major");
    test(mra, mrb, mrc, "Hybrid row-major");
    test(mrans, mcbns, mrcns, "Hybrid col-major and row-major, no shark tooth");
    test(mrans, mrbns, mrcns, "Hybrid row-major, no shark tooth");
    test(mraf, mcbf, mrcf, "Hybrid col-major and row-major with float");
    test(mra, mcb, mrc, "Hybrid col-major and row-major");
    test(mzra, mzcb, mzrc, "Hybrid col-major and row-major, Z-order");
    test(mra, mzcb, mzrc, "Hybrid col-major and row-major, Z and E-order");
    test(mra, dcb, mzrc, "Hybrid col-major and row-major, Z and E-order mixed with dense2D");
    test(mra, db, mrcns, "Hybric matrix = Shark * dense2D");
    test(mrans, db, mccns, "Hybric matrix (col-major) = hybrid (row) * dense2D");

    return 0;
}
 














