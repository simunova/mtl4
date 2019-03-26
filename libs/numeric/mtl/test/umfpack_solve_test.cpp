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

// Example build call: scons -D check=1 with_blas=1 blas_path=/usr/lib blas_ldflags=-lblas with_umfpack=1 umfpack_path=/u/pgottsch/Software/UMFPACK-5.3.0 amd_path=/u/pgottsch/Software/AMD ufconfig_path=/u/pgottsch/Software/UFconfig umfpack_solve_test

// Example build call: scons -D check=1 with_blas=1 blas_path=/usr/lib64 blas_ldflags=-lblas with_umfpack=1 umfpack_path=/home/wr2/pgottsch/64bit/projects/amdis/AMDiS/lib/UMFPACK amd_path=/home/wr2/pgottsch/64bit/projects/amdis/AMDiS/lib/AMD ufconfig_path=/home/wr2/pgottsch/64bit/projects/amdis/AMDiS/lib/UFconfig umfpack_solve_test

// CMAKE_CXX_FLAGS: -DMTL_HAS_UMFPACK -I/home/pgottsch/Software/UMFPACK-5.3.0/Include -I/home/pgottsch/Software/AMD/Include -I/home/pgottsch/Software/UFconfig

// CMAKE_EXE_LINKER_FLAGS: -L/home/pgottsch/Software/UMFPACK-5.3.0/Lib -lumfpack -L/home/pgottsch/Software/AMD/Lib -lamd -lblas (are before umfpack_solve_test.o but must come after! -> fixed link.flags)

#include <iostream>
#include <cmath>
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>
 

using namespace std;  

inline void add_imag(float&, double) {}
inline void add_imag(double&, double) {}
inline void add_imag(long double&, double) {}
inline void add_imag(complex<float>& v, double inc) { v+= complex<float>(0, inc); }
inline void add_imag(complex<double>& v, double inc) { v+= complex<double>(0, inc); }
inline void add_imag(complex<long double>& v, double inc) { v+= complex<double>(0, inc); }

#ifdef MTL_HAS_UMFPACK
template <typename Matrix>
int test(const Matrix&, const char* name)
{
    using mtl::Collection;
    
    typedef typename Collection<Matrix>::value_type   value_type;
    
    value_type array[5][5]= {{ 2.,  3.,  0.,  0.,  0.},
			     { 3.,  0.,  4.,  0.,  6.},
			     { 0., -1., -3.,  2.,  0.},
			     { 0.,  0.,  1.,  0.,  0.},
			     { 0.,  4.,  2.,  0.,  1.}};
    Matrix A(array);
    crop(A);
    if (A.nnz() != 12) 
	throw "Matrix should have 12 non-zeros!";

    value_type                      b_array[5]= {8., 45., -3., 3., 19.};
    mtl::dense_vector<value_type>   x(5), b(b_array);

    cout << name << "\nA = \n" << A << "b = " << b << "\n";

    mtl::mat::umfpack::solver<Matrix> solver(A);
    int status= solver(x, b);
    // int status= umfpack_solve(A, x, b); // creates solver on the fly
    cout << "A \\ b = " << x << "\n\n";

    for (int i= 0; i < 5; i++) 
	if (std::abs(x[i] - value_type(i+1)) > 0.01)
	    throw "Wrong result!";

    {
	mtl::mat::inserter<Matrix> ins(A);
	value_type v(5); add_imag(v, 1.); // set to 5 or 5+i depending on type
	ins[1][2] << v;
    }
    b[1]= 48.; add_imag(b[1], 3.); // set to 48 or 48+3i depending on type
    cout << "\nA = \n" << A << "b = " << b << "\n";
    solver.update_numeric();

    status= solver(x, b);
    cout << "A \\ b = " << x << "\n\n";

    for (int i= 0; i < 5; i++) 
	if (std::abs(x[i] - value_type(i+1)) > 0.01)
	    throw "Wrong result after update_numeric!";

    {
	mtl::mat::inserter<Matrix> ins(A);
	ins[3][4] << 2.;
    }
    b[3]= 13.;
    cout << "\nA = \n" << A << "b = " << b << "\n";
    solver.update();

    status= solver(x, b);
    cout << "A \\ b = " << x << "\n\n";

    for (int i= 0; i < 5; i++) 
	if (std::abs(x[i] - value_type(i+1)) > 0.01)
	    throw "Wrong result after update!";

    {
	mtl::mat::inserter<Matrix> ins(A);
	ins[3][4] << 3.;
    }
    b[3]= 18.;
    cout << "\nA = \n" << A << "b = " << b << "\n";

    // creates solver on the fly
    status= umfpack_solve(A, x, b);
    cout << "A \\ b = " << x << "\n\n";

    for (int i= 0; i < 5; i++) 
	if (std::abs(x[i] - value_type(i+1)) > 0.01)
	    throw "Wrong result after update!";
    return status;
}
#endif


int main(int, char**)
{
#ifdef MTL_HAS_UMFPACK
    using namespace mtl;
    typedef mat::parameters<col_major>           col_para;

#if 0 // weird error, will be fixed if someone really uses this
    test(compressed2D<complex<long double> >(),          "complex<long double> row-major");
    test(compressed2D<complex<long double>, col_para>(), "complex<long double> column-major");
#endif

    test(compressed2D<long double>(),          "long double row-major");
    test(compressed2D<long double, col_para>(), "long double column-major");

    test(compressed2D<complex<double> >(),          "complex<double> row-major");
    test(compressed2D<complex<double>, col_para>(), "complex<double> column-major");

    test(compressed2D<complex<float> >(),           "complex<float> row-major");
    test(compressed2D<complex<float>, col_para>(),  "complex<float> column-major");

    test(compressed2D<double>(),                    "double row-major");
    test(compressed2D<double, col_para>(),          "double column-major");

    test(compressed2D<float>(),                     "float row-major");
    test(compressed2D<float, col_para>(),           "float column-major");

#else
    std::cout << "Test is ignored when MTL_HAS_UMFPACK is not defined\n";
#endif

    return 0;
}
