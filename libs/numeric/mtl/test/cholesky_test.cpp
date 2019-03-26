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
#include <cmath>

#include <boost/numeric/mtl/mtl.hpp>


using namespace std;  
using mtl::generate_mask; using mtl::row_major; using mtl::col_major;

template <typename Matrix>
void print_matrix(Matrix& matrix)
{ 
    using std::cout;
    typedef typename mtl::Collection<Matrix>::size_type size_type;
    for (size_type i= 0 ; i < num_rows(matrix); i++ ){
	for(size_type j=0; j < num_cols(matrix);  j++ ){
	    cout.fill (' '); cout.width (8); cout.precision (5); cout.flags (ios_base::left);
	    cout << showpoint <<  matrix[i][j] <<"  ";
	}
	cout << endl;
    }
}



template <typename Matrix>
void test(Matrix& matrix, const char* name)
{
    namespace with_bracket = mtl::mat::with_bracket;
    // namespace with_iterator = mtl::mat::with_iterator;

    using mtl::mat::recursive_cholesky_visitor_t;
    using mtl::mat::detail::mult_schur_update_t;

    std::cout << "Test " << name << "\n-----\n\n";
    fill_matrix_for_cholesky(matrix);

    recursive_cholesky(matrix);
    if (matrix.num_cols() <= 10) { 
	print_matrix(matrix); std::cout << "\n"; 
    }

    fill_matrix_for_cholesky(matrix);

#if 0
    with_iterator::recursive_cholesky_base_visitor_t  iter_vis;
    recursive_cholesky(matrix, iter_vis);
    if (matrix.num_cols() <= 10) { 
	print_matrix(matrix); std::cout << "\n"; 
    }

    fill_matrix_for_cholesky(matrix);
#endif

    recursive_cholesky_visitor_t<mtl::recursion::bound_test_static<2>, with_bracket::cholesky_base_t, with_bracket::tri_solve_base_t, 
                                 with_bracket::tri_schur_base_t, with_bracket::schur_update_base_t>   
        iter_vis2; 
    recursive_cholesky(matrix, iter_vis2);
    if (matrix.num_cols() <= 10) { 
	print_matrix(matrix); std::cout << "\n"; 
    }

    fill_matrix_for_cholesky(matrix);


#if 0 // ITERATOR VERSION CURRENTLY NOT SUPPORTED -- CAUSES SEGFAULT, e.g. with icc 11.0 in r8536 !!!!

    recursive_cholesky_visitor_t<mtl::recursion::bound_test_static<2>, with_iterator::cholesky_base_t, with_iterator::tri_solve_base_t, 
                                 with_iterator::tri_schur_base_t, with_iterator::schur_update_base_t>   
        iter_vis3;

    recursive_cholesky(matrix, iter_vis3);
    if (matrix.num_cols() <= 10) { 
	print_matrix(matrix); std::cout << "\n"; 
    }


    fill_matrix_for_cholesky(matrix);

    typedef mult_schur_update_t<mtl::gen_tiling_22_dmat_dmat_mult_t<mtl::assign::minus_sum> > schur_update_22_t;
    recursive_cholesky_visitor_t<mtl::recursion::bound_test_static<2>, with_iterator::cholesky_base_t, with_iterator::tri_solve_base_t, 
                                 with_iterator::tri_schur_base_t, schur_update_22_t>   
        iter_vis4;

    recursive_cholesky(matrix, iter_vis4);
    if (matrix.num_cols() <= 10) { 
	print_matrix(matrix); std::cout << "\n"; 
    }


    typedef detail::mult_schur_update_t<gen_tiling_44_dmat_dmat_mult_t<minus_mult_assign_t> > schur_update_44_t;

#endif
}



int main(int argc, char* argv[])
{
    mtl::vampir_trace<9999> tracer;

    using namespace mtl;
    unsigned size= 13; 
    if (argc > 1) size= atoi(argv[1]); 

    dense2D<double>                                dr(size, size);
    dense2D<double, mat::parameters<col_major> > dc(size, size);
    morton_dense<double,  morton_mask>             md(size, size);
    morton_dense<double,  morton_z_mask>           mzd(size, size);
    morton_dense<double,  doppled_2_row_mask>      d2r(size, size);
    morton_dense<double,  doppled_2_col_mask>      d2c(size, size);
    morton_dense<double,  doppled_16_row_mask>     d16r(size, size);
    morton_dense<double,  doppled_32_row_mask>     d32r(size, size);
    morton_dense<double,  doppled_64_row_mask>     d64r(size, size);
    morton_dense<double,  doppled_64_col_mask>     d64c(size, size);
    morton_dense<double,  doppled_128_col_mask>    d128r(size, size);
    size= 9; 
    dense2D<double>                                dr2(size, size);
    
    test(dr2, "Dense row major");
    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(md, "Morton N-order");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(d2c, "Hybrid 2 column-major");
    test(d16r, "Hybrid 16 row-major");


    return 0;
}





