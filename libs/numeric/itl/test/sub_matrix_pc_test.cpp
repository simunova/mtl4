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

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

template <typename Matrix>
inline void strided_laplacian_setup(Matrix& A, unsigned m, unsigned n)
{
    A.change_dim(2*m*n, 2*m*n);
    set_to_zero(A);
    mtl::mat::inserter<Matrix>      ins(A, 5);

    for (unsigned i= 0; i < m; i++)
	for (unsigned j= 0; j < n; j++) {
	    typename mtl::Collection<Matrix>::value_type four(4.0), minus_one(-1.0);
	    unsigned row= 2 * (i * n + j);
	    ins(row, row) << four;
	    if (j < n-1) ins(row, row+2) << minus_one;
	    if (i < m-1) ins(row, row+2*n) << minus_one;
	    if (j > 0) ins(row, row-2) << minus_one;
	    if (i > 0) ins(row, row-2*n) << minus_one;
	}
    for (unsigned i= 1; i < num_rows(A); i+= 2)
	ins(i, i) << 2;
} 


int main()
{
    using mtl::srange; using mtl::imax;
    // For a more realistic example set sz to 1000 or larger
    const int size = 3, N = 2 * size * size; 

    typedef mtl::compressed2D<double>  matrix_type;
    typedef itl::pc::ic_0<matrix_type> ic_type;


    mtl::compressed2D<double>          A;
    strided_laplacian_setup(A, size, size);
    mtl::io::tout << "A (merged diagonal and Laplace) is\n" << A << '\n';

    mtl::dense_vector<bool> tags= make_tag_vector(N, srange(0, imax, 2));
    itl::pc::sub_matrix_pc<ic_type, matrix_type> P(tags, A);

    mtl::dense_vector<double>          x(N, 1.0), b(N);
    
    b = A * x;
    x= 0;

    itl::cyclic_iteration<double> iter(b, N, 1.e-6, 0.0, 3);
    cg(A, x, b, P, iter);
     
    // Test if adjoint works
    x= 0;
    itl::cyclic_iteration<double> iter2(b, N, 1.e-6, 0.0, 3);
    bicg(A, x, b, P, iter2);

    return 0;
}
