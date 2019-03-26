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

// With contributions from Cornelius Steinhardt

#include <cstdlib>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/operation/cuppen.hpp>

using namespace std;

const double 			tol= 1.0e-5;

template <typename Matrix, typename Value, typename Vector>
void test_vector(const Matrix& A, const Value& alpha, const Vector& v, int i)
{
    Vector diff(A*v-alpha*v), v1(A*v), v2(alpha*v); // , diff(v1-v2);
    if (size(v1) < 17) 	
	cout << "A*v is     " << v1 << "\nalpha*v is " << v2 << '\n';
    if (two_norm(diff) > tol) cout << "two_norm(difference) of the " << i << "-th eigenvector is " << two_norm(diff) << '\n'; // throw "wrong eigenvector";
}

template <typename Matrix, typename Value, typename Vector>
void test(const Matrix& B, Matrix& BQ, Value scaling, Vector& lambda_b)
{
    if (num_rows(B) <= 20)
	std::cout << "B=\n" << B << "\n";

    mtl::dense_vector<double> eig_b= eigenvalue_symmetric(B,22);
    sort(eig_b);
    eig_b*= scaling;
    std::cout<<"eigenvalues with QR ="<< eig_b <<"\n";
    
    cuppen(B, BQ, lambda_b);
    lambda_b*= scaling;
    // std::cout<<"B  =\n"<< B <<"\n";
    if (num_rows(B) <= 20)
	std::cout<<"Q  =\n"<< BQ <<"\n";
    std::cout<<"eigenvalues with Cuppen ="<< lambda_b <<"\n";
    
    eig_b-= lambda_b;
    std::cout<<"two_norm(diff)  ="<< two_norm(eig_b) <<"\n";
    MTL_THROW_IF(two_norm(eig_b) > tol, mtl::runtime_error("Cuppen computes wrong eigenvalues"));

    Matrix Bs(scaling * B);
    for (unsigned i= 0; i < num_rows(B); i++)
	test_vector(Bs, lambda_b[i], mtl::dense_vector<double>(BQ[mtl::iall][i]), i);
}

int main(int argc, char** argv)
{
    using namespace mtl;
    int size= 16;
    if (argc > 1) size= atoi(argv[1]);

    dense_vector<double>        eig, lam(2), lambda(4),  lambda_b(size), eig_b(size);
#if 1
    dense2D<double> Mini(2, 2), QM(2, 2);
    Mini= 1, 2, 2, -11;
    cuppen(Mini, QM, lam);
    
    std::cout << "eigenvalues of Mini  =" << lam <<"\n";
    std::cout << "eigenvectors of Mini are\n" << QM;    
#endif

    double array[][4]= {{1,  2,   0,  0},
                        {2, -9,  -2,  0},
                        {0, -2,   1,  3},
                        {0,  0,   3, 10}};
    dense2D<double> A(array), Q(4,4);
    std::cout << "A=\n" << A << "\n";

    eig= eigenvalue_symmetric(A,22);
    sort(eig);
    std::cout<<"eigenvalues  ="<< eig <<"\n";
    
    cuppen(A, Q, lambda);
    std::cout<<"A  =\n"<< A <<"\n";
    std::cout<<"Q  =\n"<< Q <<"\n";
    std::cout<<"eigenvalues  ="<< lambda <<"\n";
   

    eig-= lambda;
    std::cout<<"two_norm(diff)  ="<< two_norm(eig) <<"\n";
    MTL_THROW_IF(two_norm(eig) > tol, mtl::runtime_error("Cuppen computes wrong eigenvalues"));

    for (unsigned i= 0; i < num_rows(A); i++)
	test_vector(A, lambda[i], dense_vector<double>(Q[iall][i]), i);

    
    dense2D<double> B(size,size), BQ(size,size);
    B= 0; BQ= 0;
    
    const double scale= 64.0 / double(size);
    for(int i= 1; i < size ; i++){
      B[i][i]= scale*i+6;
      B[i][i-1]= 1;
      B[i-1][i]= 1;
    }
    B[0][0]= 6;

    test(B, BQ, 1.0, lambda_b);
    
    const double maxv= B[size-1][size-1];
    B/= maxv;
    test(B, BQ, maxv, lambda_b);

#if 0
    // Poisson equation cannot be solved, double eigenvalues are now correctly handled by the secular equation
    // but Q_tilde in cuppen contains nans (0/0)
    int lsize= 4;
    if (argc > 1) lsize= atoi(argv[1]);

    dense2D<double> C(lsize, lsize), CQ(lsize, lsize);
    C= 0; CQ= 0;
    
    C[0][0]= 2;
    for(int i= 1; i < lsize; i++) {
	C[i][i]= 2;
	C[i][i-1]= -1;
	C[i-1][i]= -1;
    }
    cout << "The matrix of the 1D-Poisson equations I\n" << C << '\n';
	

    dense_vector<double> lambda_c(lsize);
    cuppen(C, CQ, lambda_c);

    if (lsize <= 100)
	cout << "The eigenvalues of the 1D-Poisson equations are " << lambda_c << '\n';
    if (lsize <= 20)
	cout << "The eigenvectors of the 1D-Poisson equations are\n" << CQ << '\n';

    for (unsigned i= 0; i < num_rows(C); i++);
	//test_vector(C, lambda_c[i], dense_vector<double>(CQ[iall][i]));
#endif
    
    return 0;
}



