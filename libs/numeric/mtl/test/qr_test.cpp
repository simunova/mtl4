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
#include <boost/utility.hpp>
#include <boost/numeric/mtl/mtl.hpp>


using namespace std;


double f(double) { /* cout << "double\n"; */ return 1.0; }
complex<double> f(complex<double>) 
{ 
    //cout << "complex\n"; 
    return complex<double>(1.0, -1.0); 
}

 
int main(int, char**)
{
    using namespace mtl; using mtl::io::tout;
    unsigned size=4, row= size+1, col=size;

    double tol(0.00001);
    dense_vector<double>                vec(size), vec1(size);
    dense2D<double>                     A(row, col),   Q(row, row),   R(row, col),   A_test(row, col),
	                                A_t(col, row), Q_t(col, col), R_t(col, row), A_t_test(col, row);
    dense2D<complex<double> >                            dz(row, col), Qz(row, row), Rz(row, col);
    dense2D<double, mat::parameters<col_major> >         dc(size, size);
    compressed2D<double>                                 Ac(size, size), Qc(size, size), Rc(size, size), A_testc(size, size) ;
    A= 0; 

    A[0][0]=1;    A[0][1]=1;    A[0][2]=1;
    A[1][0]=3;    A[1][1]=-1;   A[1][2]=-2;
    A[2][0]=1;    A[2][1]=7;    A[2][2]=1;
    A[3][3]=-10;  A[4][0]=4;    A[4][2]=3;
    tout<<"A=\n"<< A <<"\n";
    laplacian_setup(Ac, 2,2);

    
    tout<<"START-----dense2d---------row > col\n";
  
    dense2D<double> A1(A[iall][iall]), A2(A);
    boost::tie(Q, R)= qr(A1);
    tout<<"R=\n"<< R <<"\n";
    tout<<"Q=\n"<< Q <<"\n";
    A_test= Q*R-A2;
    tout<<"Q*R=\n"<< Q*R <<"\n";
	
    tout<< "one_norm(Rest A)=" << one_norm(A_test) << "\n";
    MTL_THROW_IF(one_norm(A_test) > tol, mtl::logic_error("wrong QR decomposition of matrix A"));

	
    tout<<"START------dense2d-------row < col\n";

    A_t= trans(A);
    boost::tie(Q_t, R_t)= qr(A_t);
    tout<<"R_t=\n"<< R_t <<"\n";
    tout<<"Q_t=\n"<< Q_t <<"\n";
    A_t_test= Q_t*R_t-A_t;
    tout<<"Q_t*R_t=\n"<< Q_t*R_t <<"\n";
			
    tout<< "one_norm(Rest A')=" << one_norm(A_t_test) << "\n";
    MTL_THROW_IF(one_norm(A_t_test) > tol, mtl::logic_error("wrong QR decomposition of matrix trans(A)"));
	
    tout<<"START-------compressed2d-------row > col\n";
#if 1
    boost::tie(Qc, Rc)= qr(Ac);
    tout<<"R=\n"<< Rc <<"\n";
    tout<<"Q=\n"<< Qc <<"\n";
    A_testc= Qc*Rc-Ac;
    tout<<"Q*R=\n"<< Qc*Rc <<"\n";
    tout<<"A=\n"<< Ac <<"\n";
	
    tout<< "one_norm(Rest A)=" << one_norm(A_testc) << "\n";
    MTL_THROW_IF(one_norm(A_testc) > tol, mtl::logic_error("wrong QR decomposition of matrix A")); 
#endif

#if 0
    dz[0][0]=complex<double>(1.0, 0.0);
    dz[0][1]=complex<double>(1.0, 0.0);
    dz[0][2]=complex<double>(1,0);
    dz[1][0]=complex<double>(1,0);
    dz[1][1]=complex<double>(-1,0);
    dz[1][2]=complex<double>(-2,0);
    dz[2][0]=complex<double>(1,0);
    dz[2][1]=complex<double>(-2,0);
    dz[2][2]=complex<double>(1,0);
    dz[3][3]=complex<double>(-10,0);
    tout<<"MAtrix complex=\n"<< dz <<"\n";

    tout<<"START-----complex---------"<< dz[0][0] << "\n";
    //
    boost::tie(Qz, Rz)= qr(dz);
    // Rz= qr_zerl(dz).second;
    // tout<<"MAtrix  R="<< Rz <<"\n";
    // tout<<"MAtrix  Q="<< Qz <<"\n";
    // tout<<"MAtrix  A=Q*R--outside"<< Qz*Rz <<"\n";
#endif
    return 0;
}

