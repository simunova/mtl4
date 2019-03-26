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

inline void error(const char* message)
{
    std::cerr << message << '\n';
//     throw message;
}

template <typename Matrix>
void test(Matrix& A, const char* name)
{
    using mtl::Collection; using mtl::operations::update_plus;
    using mtl::mat::inserter; using mtl::mat::shifted_inserter; 

    typedef typename Collection<Matrix>::value_type   value_type;
    
    value_type array[][4]= {{3, 7.2, 0, 4}, {2, 4.444, 5, 2}, {1, 3, 8, 1}};
    A= array;

    std::cout << "\n" << name << ", assignment: A = \n" << A << "\n";

    {
	shifted_inserter< inserter<Matrix, update_plus<value_type> > > ins(A, 0, 1, 2);

	ins[0][0] << 3;
	ins[1][0]+= -1;
	ins[1][1]= 1.5;
    }
    std::cout <<  "... after shifted insertion A = \n" << A << "\n";
    
    if (num_rows(A) != 3 || num_cols(A) != 4) error("Wrong matrix size");

    if (A[1][0] != value_type(2)) error("Value should be unchanged");

    if (A[1][2] != value_type(8)) error("Error in shifted updating");

    if (A[2][2] != value_type(7)) error("Error in shifted +=");

    if (A[2][3] != value_type(1.5)) error("Error in shifted =");

    if (A[1][3] != value_type(2)) error("Value should be unchanged");
    
    if (A[2][2] != value_type(8)) error("Error ok. Test only for coverage"); 
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
