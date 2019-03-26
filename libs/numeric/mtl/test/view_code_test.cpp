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
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/view_code.hpp>

using namespace std;
typedef std::complex<double> ct;

// This template should never be instantiated
template <int Code, typename Matrix>
struct check_type_aux
{
    static void apply() { MTL_STATIC_ASSERT((Code == 999), "View not folded correctly."); }
};

template <typename Value, typename Parameters>
struct check_type_aux<0, mtl::mat::dense2D<Value, Parameters> >
{    static void apply() {}  };

template <typename Matrix>
struct check_type_aux<2, mtl::mat::conj_view<Matrix> > 
{    static void apply() {}  };

template <typename Matrix>
struct check_type_aux<4, mtl::mat::transposed_view<Matrix> > 
{    static void apply() {}  };

template <typename Matrix>
struct check_type_aux<6, mtl::mat::hermitian_view<Matrix> > 
{    static void apply() {}  };


template <int Code, typename Matrix>
void check_code()
{
    // cout << "View code is " << mtl::traits::view_code<Matrix>::value << '\n';
    // MTL_THROW_IF((Code != (mtl::traits::view_code<Matrix>::value & 6)), mtl::unexpected_result("Not correct view_code."));
    MTL_STATIC_ASSERT((Code == (mtl::traits::view_code<Matrix>::value & 6)), "Not correct view_code.");
}

template <int Code, typename Matrix>
void check_type(const Matrix&)
{
    check_code<Code, Matrix>();
    check_type_aux<Code, Matrix>::apply();
}


template <typename Matrix>
void check(const Matrix& A, int code)
{
    switch (code) {
      case 0: MTL_THROW_IF(A[1][0] != ct(3, 4), mtl::unexpected_result("In original matrix.")); break;
      case 2: MTL_THROW_IF(A[1][0] != ct(3, -4), mtl::unexpected_result("In conjugated matrix.")); break;
      case 4: MTL_THROW_IF(A[1][0] != ct(0, 1), mtl::unexpected_result("In transposed matrix.")); break;
      case 6: MTL_THROW_IF(A[1][0] != ct(0, -1), mtl::unexpected_result("In hermitian matrix.")); break;
    }
}

template <int Code, typename Matrix>
void inner_test(const Matrix& A, const char* text)
{
    cout << "\nChecking " << text << ", A is\n" << A;
    check(A, Code);
    check_type<Code>(A);

    cout << "conj(A) is\n" << mtl::conj(A);
    check(mtl::conj(A), Code ^ 2);
    check_type<Code ^ 2>(mtl::conj(A));

    cout << "trans(A) is\n" << trans(A);
    check(trans(A), Code ^ 4);
    check_type<Code ^ 4>(trans(A));

    cout << "hermitian(A) is\n" << hermitian(A);
    check(hermitian(A), Code ^ 6);
    check_type<Code ^ 6>(hermitian(A));
}


template <typename Matrix>
void test(Matrix& A)
{
    A= ct(2, 1), ct(0, 1),
       ct(3, 4), ct(5, 6);

    inner_test<0>(A, "original matrix");
    inner_test<2>(mtl::conj(A), "conjugated matrix");
    inner_test<4>(trans(A), "transposed matrix");
    inner_test<6>(hermitian(A), "Hermitian matrix");
}



int main(int, char**)
{
    using namespace mtl;

    dense2D<ct>    A(2, 2);

    test(A);

    return 0;
}
