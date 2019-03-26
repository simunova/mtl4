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

#define MTL_WITH_DEVELOPMENT
#define MTL_VERBOSE_TEST

#include <string>
#include <iostream>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/matrix/sparse_banded.hpp>


template <typename Matrix>
void laplacian_test(Matrix& A, unsigned dim1, unsigned dim2, const char* name)
{
    mtl::io::tout << "\n" << name << "\n";
    laplacian_setup(A, dim1, dim2);
    mtl::io::tout << "Laplacian A:\n" << A << std::endl;
    if (dim1 > 1 && dim2 > 1) {
	typename Matrix::value_type four(4.0), minus_one(-1.0), zero(0.0);
	MTL_THROW_IF(A[0][0] != four, mtl::runtime_error("wrong diagonal"));
	MTL_THROW_IF(A[0][1] != minus_one, mtl::runtime_error("wrong east neighbor"));
	MTL_THROW_IF(A[0][dim2] != minus_one, mtl::runtime_error("wrong south neighbor"));
	MTL_THROW_IF(dim2 > 2 && A[0][2] != zero, mtl::runtime_error("wrong zero-element"));
	MTL_THROW_IF(A[1][0] != minus_one, mtl::runtime_error("wrong west neighbor"));
	MTL_THROW_IF(A[dim2][0] != minus_one, mtl::runtime_error("wrong north neighbor"));
	MTL_THROW_IF(dim2 > 2 && A[2][0] != zero, mtl::runtime_error("wrong zero-element"));
    }
}

template <typename Matrix>
void rectangle_test(Matrix& A, const char* name)
{
    {
	mtl::mat::inserter<Matrix> ins(A);
	int i= 1;
	unsigned nc= num_cols(A);
	for (unsigned r= 0; r < num_rows(A); r++) {
	    if (r < nc - 4) ins(r, r + 4) << i++;
	    if (r < nc) ins(r, r) << i++;
	    if (r >= 2 && r < nc + 2) ins(r, r - 2) << i++;
	    if (r >= 4 && r < nc + 4) ins[r][r - 4] << i++;
	}
    }
    mtl::io::tout << name << ": A=\n" << A << '\n';
}

template <typename Matrix, typename Tag>
void two_d_iteration(const Matrix & A, Tag)
{
    namespace traits = mtl::traits;

    typename traits::row<Matrix>::type                                 row(A); 
    typename traits::col<Matrix>::type                                 col(A); 
    typename traits::const_value<Matrix>::type                         value(A); 
    typedef typename traits::range_generator<Tag, Matrix>::type        cursor_type;
    for (cursor_type cursor = mtl::begin<Tag>(A), cend = mtl::end<Tag>(A); cursor != cend; ++cursor) {
	typedef mtl::tag::nz     inner_tag;
	mtl::io::tout << "---\n";
	typedef typename traits::range_generator<inner_tag, cursor_type>::type icursor_type;
	for (icursor_type icursor = mtl::begin<inner_tag>(cursor), icend = mtl::end<inner_tag>(cursor); icursor != icend; ++icursor)
	    mtl::io::tout << "A[" << row(*icursor) << ", " << col(*icursor) << "] = " << value(*icursor) << '\n';
    }
    mtl::io::tout << "===\n\n";
} 

template <typename Matrix>
void mat_vec_mult_test(const Matrix& A, const char* name)
{
    typedef typename Matrix::value_type  value_type;
    mtl::io::tout << name << " " << num_rows(A) << " by " << num_cols(A) << '\n' << A;

    mtl::dense_vector<value_type> v, w(num_cols(A), 3.0), v2;
    v= A * w;
    mtl::io::tout << "A * v =    " << v << '\n';

    mtl::compressed2D<value_type> B(A);
    v2= B * w;
    mtl::io::tout << "Should be: " << v2 << "\n\n";
    v2-= v;
    MTL_THROW_IF(two_norm(v2) > 0.001, 
		 mtl::runtime_error("wrong result for sparse banded times vector"));
}

int main(int, char**)
{
    using namespace mtl;
#ifdef MTL_WITH_DEVELOPMENT
    unsigned dim1= 3, dim2= 4;
    mat::sparse_banded<double>  dr, dr2(6, 11), dr3(11, 6), dr4(6, 5);
    
    laplacian_test(dr, dim1, dim2, "Dense row major");
    rectangle_test(dr2, "Dense row major");
    rectangle_test(dr3, "Dense row major");
    rectangle_test(dr4, "Dense row major");

    mat::compressed2D<double> C;
    laplacian_setup(C, dim1, dim2);

    mat::sparse_banded<double>  D;
    D= C;
    mtl::io::tout << "D is\n" << D << '\n';

    two_d_iteration(D, mtl::tag::row());

    mat::compressed2D<double> E;
    E= D;
    mtl::io::tout << "E is\n" << E << '\n';

    mat::sparse_banded<double>  dr5(5, 5), dr6(5, 5);
    {
	mtl::mat::inserter<mat::sparse_banded<double> > ins5(dr5), ins6(dr6);	
	ins5[2][0] << 1; ins5[3][1] << 2; ins5[4][2] << 3; ins5[4][0] << 4;
	ins6[0][2] << 1; ins6[1][3] << 2; ins6[2][4] << 3; ins6[0][4] << 4;
    }

    mat_vec_mult_test(dr2, "Dense row major");
    mat_vec_mult_test(dr3, "Dense row major");
    mat_vec_mult_test(dr4, "Dense row major");
    mat_vec_mult_test(dr5, "Dense row major");
    mat_vec_mult_test(dr6, "Dense row major");

    mat_vec_mult_test(dr, "Dense row major");
#endif
    


    return 0;
}
 
