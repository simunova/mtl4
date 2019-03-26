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
#include <cstdlib>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

void test_dynam_vector()
{
    std::cout<<"\n---test: operations for vectors of dynam. size---";
	
    typedef mtl::dense_vector<double> vector;
    vector v1(10,5.), v2(10,0.), v3(10);
	
    v1[7]=3.;
	
    v2[2]=2.;
    std::cout<<"\nv1="<<v1;
    std::cout<<"\nv2="<<v2;

    v3= v1 + 5. * v2 - dot(v1,v2) * v2;
    std::cout << "\nv3=v1+5.*v2-mtl::dot(v1,v2)*v2=" << v3;
    v3/= 5.;
    std::cout << "\nv3/=5. is " << v3;
    std::cout << "\nv3 has size=" << size(v3);
}

void test_stat_vector()
{
    std::cout<<"\n---test: operations for vectors of fixed size---";
    typedef mtl::vec::parameters<mtl::tag::col_major, mtl::vec::fixed::dimension<10> > dimension;
    typedef mtl::dense_vector<double, dimension> vector;
    vector v1, v2, v3;
    v1=5.;
    v1[7]=3.;
    v2=0.;
    v2[2]=2.;
    std::cout<<"\nv1="<<v1;
    std::cout<<"\nv2="<<v2;
    v3=v1+5.*v2-mtl::dot(v1,v2)*v2;
    std::cout<<"\nv3=v1+5.*v2-mtl::dot(v1,v2)*v2="<<v3;
    v3/=5.;
    std::cout<<"\nv3/=5. is "<<v3;
    std::cout << "\nv3 has size=" << mtl::static_size<vector>::value;
}

void test_stat_matrix()
{
    std::cout<<"\n---test: operations for matrices of stat. size---";
    typedef mtl::mat::parameters<mtl::tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3> > matrix_parameters;
    typedef mtl::dense2D<double, matrix_parameters>  matrix_type;
	
    matrix_type A, B, C;
	
    mtl::mat::diagonal_setup(A, 2. );
	
    B=5.;
    mtl::mat::inserter<matrix_type> ins(B);
    ins[1][0]<<2.;
	
    B(0,1)=3.;
	
    std::cout<<"\nA="<<A;
    std::cout<<"\nB="<<B;
    std::cout<<"\nmatrices A,B,C have size=("<<mtl::static_num_rows<matrix_type>::value<<","<<mtl::static_num_cols<matrix_type>::value<<")";
    std::cout<<"\nmatrix B has one-norm ="<<mtl::one_norm(B);
	
    invert_diagonal(A);
    std::cout<<"\ninverto of diagonal matrix A is="<<A;

    C=A+B;	
    B += A;
    C = 5.*B -A;
    C=A*B;
    mtl::dense_vector<mtl::Collection<matrix_type>::value_type> eigenvalues;
    eigenvalues= eigenvalue_symmetric(A);
    std::cout<<"\nmatrix A has trace ="<<trace(A);
}



void test_dynam_matrix()
{
    std::cout<<"\n---test: operations for matrices of dynam. size---";
    typedef mtl::dense2D<double>  MATRIX;
	
    MATRIX A(3,3),B(3,3),C(3,3);
	
    mtl::mat::diagonal_setup(A, 2. );
	
    B=5.;
    mtl::mat::inserter<MATRIX> ins(B);
    ins[1][0]<<2.;
	
    B(0,1)=3.;
	
    std::cout<<"\nA="<<A;
    std::cout<<"\nB="<<B;
    std::cout<<"\nMatrix A has size=("<<A.num_rows()<<","<<A.num_cols()<<")";
    std::cout<<"\nmatrix B has one-norm ="<<mtl::one_norm(B);
	
    invert_diagonal(A);
    std::cout<<"\ninvert of diagonal matrix A is="<<A;

    C=A+B;
    B += A; 
    C = 5.*B -A;
    C=A*B;                         
    mtl::dense_vector<mtl::Collection<MATRIX>::value_type > eigenvalues;
    eigenvalues=eigenvalue_symmetric(A);
}


void test_compressed_matrix()
{

    std::cout<<"\n---test: operations for compressed matrices ---";
    typedef mtl::compressed2D<double>  MATRIX;
	
    MATRIX A(3,3),B(3,3),C(3,3);
	
    mtl::mat::diagonal_setup(A, 2. );
	
    {
	mtl::mat::inserter<MATRIX> ins(B);
	ins[1][0]<<2.;
	ins[0][1]<<3.;
    }
	
    std::cout<<"\nA="<<A;
    std::cout<<"\nB="<<B;
    std::cout<<"\nMatrix A has size=("<<A.num_rows()<<","<<A.num_cols()<<")";
    std::cout<<"\nmatrix B has one-norm ="<<mtl::one_norm(B);
	
    invert_diagonal(A);
    std::cout<<"\ninverto of diagonal matrix A is="<<A;
	
    C=A+B;
    std::cout<<"\nC=A+B="<<C;
	 
    B += A;
    std::cout<<"\nB+=A; B="<<B;
	 
    C = 5.*B -A;
    std::cout<<"\nC = 5.*B -A="<<C;
	 
    C=A*B;
    std::cout<<"\nC = A*B="<<C;
	 
    mtl::dense_vector<mtl::Collection<MATRIX>::value_type > eigenvalues;
    eigenvalues= eigenvalue_symmetric(A);
}


void test_dynam_vector_and_matrix()
{

    std::cout<<"\n---test: operations for vectors and matrices of dynam. size---";
	
    typedef mtl::dense_vector<double> vector;
    vector v1(10,5.), v2(10,0.), v3(10);
    v1[7]=3.;
	
    typedef mtl::dense2D<double> matrix;
    matrix A(10,10);
    A=10.;

    v2+=A*v1;
    v3=A*v2;
}



void test_stat_vector_and_matrix()
{
    std::cout<<"\n---test: operations for vectors and matrices of fixed size---";
	
    typedef mtl::vec::parameters<mtl::tag::col_major, mtl::vec::fixed::dimension<3> > dimension;
    // typedef mtl::parameters<mtl::tag::col_major, mtl::fixed::dimension<10> > dimension;
    typedef mtl::vec::dense_vector<double, dimension> vector;
    vector v1,v2,v3;
    v1=2.; v2=0.;
		
    typedef mtl::mat::parameters<mtl::tag::row_major, mtl::index::c_index, mtl::fixed::dimensions<3, 3> > matrix_parameters;
    typedef mtl::dense2D<double, matrix_parameters>  matrix;

    matrix A;
    A=10.;
    v2+=A*v1;                            
    v3=A*v2;                    
}


void test_compressed_matrix_and_vector()
{

    std::cout<<"\n---test: operations for vectors and compressed matrices ---";
	
    typedef mtl::compressed2D<double>  matrix;
    matrix A(10,10);
	
    mtl::mat::diagonal_setup(A, 10. );
	
    typedef mtl::vec::parameters<mtl::tag::col_major, mtl::vec::fixed::dimension<10> > dimension;
    typedef mtl::vec::dense_vector<double, dimension> vector;
    vector v1,v2,v3;
    v1=2.;v2=0.;
    std::cout<<"\nA="<<A;
    std::cout<<"\nv1="<<v1;
    std::cout<<"\nv2="<<v2;
    v2+=A*v1;
    std::cout<<"\nv2+=A*v1; v2="<<v2;      
	
    v3=A*v2;
    std::cout<<"\nv3=A*v2="<<v3;
}

void test_solver()
{
    std::cout<<"\n---test: solver CG with iLU on compresed matrix ---";
    using namespace mtl;
    using namespace itl;
    const int size=5, N=size*size;

	
    typedef compressed2D<double> MATRIX;
    MATRIX A(N,N);
    mat::laplacian_setup(A,size,size);
	
    typedef dense_vector<double> vector; 
    vector x(N,1.),b(N);
    b=A*x; x=0.;
    std::cout<<"\nA="<<A;
    std::cout<<"\nb="<<b;

    pc::ilu_0<MATRIX>  Precond(A);
    noisy_iteration<double> iter(b, 500, 1.e-6);
    cg(A, x, b, Precond, iter);
    if(iter.error_code()!=0) {
	std::cout<<"\nInterpolation matrix:\n"<<A;
	std::cerr<<"\nERROR: unsolvable system, error code="<<iter.error_code(); exit(EXIT_FAILURE);
    }
    std::cout<<"\nFor compressed matrix A="<<A;
    std::cout<<"\nCG solver with iLU, b="<<b;
}

int main()
{
    test_dynam_vector();
    test_dynam_matrix();
    test_stat_vector();
    test_stat_matrix();
    test_dynam_matrix();
    test_compressed_matrix();
    test_dynam_vector_and_matrix();
    test_stat_vector_and_matrix();
    test_compressed_matrix_and_vector();
    test_solver();
    std::cout << "\nNo errors detected.\n";

    return 0;
}
