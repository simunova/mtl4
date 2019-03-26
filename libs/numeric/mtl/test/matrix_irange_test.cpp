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
#include <boost/numeric/mtl/recursion/matrix_recursator.hpp>


using namespace std;  


template <typename Size, typename Derived>
struct base1
{
    int operator[](Size) { cout << "Size\n"; return 1; }
    int operator[](mtl::irange i) 
    { 
	cout << "irange\n"; 
	return static_cast<Derived*>(this)
	    ->susu(i.start(), i.finish());
    }
};
    

template <typename Size>
struct ss : public base1<Size, ss<Size> >
{
    int susu(int, int) { cout << "susu\n"; return 3; }
};

template <typename Size>
struct ss2 : public base1<Size, ss<Size> >
{};


// For Morton matrices not applicable
template <typename Matrix>
void test(Matrix& A, const char* name)
{
    using mtl::irange; using mtl::imax;

    A= 0.0;
    A[1][1]= 1.0; 
    hessian_setup(A, 1.0);

    std::cout << "\n" << name << "\nA == \n" << A;
    
    cout << "A[irange(1, 4)][irange(1, imax)] == \n" 
	 << A[irange(1, 4)][irange(1, imax)] << "\n";

    Matrix B(A[irange(1, 4)][irange(1, imax)]);
    MTL_THROW_IF(B[1][1] != 4.0, mtl::runtime_error("Wrong value in B"));

    MTL_THROW_IF(A[irange(1, 4)][irange(1, imax)][1][1] != 4.0, mtl::runtime_error("Wrong value in A[][]"));

    Matrix C(A[irange(3, imax)][irange(1, 2)]);
    std::cout << "\n" << name << "\nA[irange(3, imax)][irange(1, 2)] == \n" << C;
    MTL_THROW_IF(C[1][0] != 5.0, mtl::runtime_error("Wrong value in C"));

    cout << "A[irange(1, 4)][irange(1, imax)][irange(1, imax)][irange(1, imax)][0][0] == \n" 
	 << A[irange(1, 4)][irange(1, imax)][irange(1, imax)][irange(1, imax)][0][0] << "\n"; 
}

template <typename Matrix>
void test2(Matrix& A, const char* name)
{
    using mtl::irange; using mtl::imax; using mtl::iall;

    A= 0.0;
    A[1][1]= 1.0; 
    hessian_setup(A, 1.0);

    std::cout << "\n" << name << "\nA == \n" << A;
    
    cout << "A[irange(2, 4)][irange(2, imax)] == \n" 
	 << A[irange(2, 4)][irange(2, imax)] << "\n";

    Matrix B(A[irange(2, 4)][irange(2, imax)]);
    MTL_THROW_IF(B[0][0] != 4.0, mtl::runtime_error("Wrong value in B"));

    MTL_THROW_IF(A[irange(2, 4)][irange(2, imax)][0][0] != 4.0, mtl::runtime_error("Wrong value in A[][]"));

    Matrix C(A[irange(4, imax)][irange(0, imax)]);
    std::cout << "\n" << name << "\nA[irange(4, imax)][irange(0, imax)] == \n" << C;
    MTL_THROW_IF(C[0][1] != 5.0, mtl::runtime_error("Wrong value in C"));

    Matrix D(A[irange(4, imax)][iall]);
    std::cout << "\n" << name << "\nA[irange(4, imax)][iall] == \n" << C;
    MTL_THROW_IF(D[0][1] != 5.0, mtl::runtime_error("Wrong value in D"));

    
}


int main(int, char**)
{
    using namespace mtl;
    const unsigned size= 5; 

    dense2D<double> dc(size, size-2);
    dense2D<double, mat::parameters<col_major> >  dcc(size, size-2);
    dense2D<float>                                   fc(size, size-2);
    morton_dense<double,  morton_mask>               mdc(size, size-2);
    morton_dense<double, doppled_32_col_mask>        mcc(size, size-2);

    test(dc, "dense2D");
    test(dcc, "dense2D col-major");

    test2(dc, "dense2D");
    test2(dcc, "dense2D col-major");
    test2(mdc, "pure Morton");
    test2(mcc, "Hybrid col-major");

    return 0;
}
