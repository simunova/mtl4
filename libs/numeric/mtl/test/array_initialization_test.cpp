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


using namespace std;


template <typename Matrix>
void test(Matrix& A, const char* name)
{
    typedef typename mtl::Collection<Matrix>::value_type   value_type;
    
#if 0 // for debuggin type infos
    typedef value_type a_type[2][3];
    typedef typename mtl::ashape::ashape<a_type>::type s_type;
    std::cout << typeid(s_type).name() << "\n";

    std::cout << boost::is_same<s_type, mtl::ashape::scal>::value  << "\n";
#endif

    value_type array[][3]= {{3, 7.2, 0}, {2, 4.444, 5}};
    A= array;

    std::cout << "\n" << name << ", assignment: A = \n" << A << "\n";

    MTL_THROW_IF(num_rows(A) != 2 || num_cols(A) != 3, mtl::runtime_error("Wrong matrix size"));
    MTL_THROW_IF(A[1][0] != value_type(2), mtl::runtime_error("Wrong value inserted"));

    Matrix B(array);

    std::cout << "\n" << name << ", construction: B = \n" << B << "\n";

    MTL_THROW_IF(num_rows(B) != 2 || num_cols(B) != 3, mtl::runtime_error("Wrong matrix size"));
    MTL_THROW_IF(B[1][0] != value_type(2), mtl::runtime_error("Wrong value inserted"));
}


int main(int, char**)
{
    using namespace mtl;
    dense2D<double>                                      dr;
    dense2D<double, mat::parameters<col_major> >      dc;
    morton_dense<double, recursion::morton_z_mask>       mzd;
    morton_dense<double, recursion::doppled_2_row_mask>  d2r;
    compressed2D<double>                                 cr;
    compressed2D<double, mat::parameters<col_major> > cc;

    dense2D<complex<double> >                            drc;
    compressed2D<complex<double> >                       crc;

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(mzd, "Morton Z-order");
    test(d2r, "Hybrid 2 row-major");
    test(cr, "Compressed row major");
    test(cc, "Compressed column major");
    test(drc, "Dense row major complex");
    test(crc, "Compressed row major complex");

	
    return 0;
}
