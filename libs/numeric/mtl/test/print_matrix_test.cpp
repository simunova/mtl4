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
#include <boost/numeric/mtl/recursion/bit_masking.hpp>
#include <boost/numeric/mtl/recursion/predefined_masks.hpp>
#include <boost/numeric/mtl/operation/print.hpp>

using namespace std;  

template <typename Matrix>
void print_matrix(Matrix& matrix)
{ 
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    using std::cout;
    for (size_type i=0 ; i<num_rows(matrix); i++ ){
	for(size_type j=0; j<num_cols(matrix);  j++ ){
	    cout.fill (' '); cout.width (8); cout.precision (5); cout.flags (ios_base::left);
	    cout << showpoint <<  matrix[i][j] <<"  ";
	}
	cout << endl;
    }
}


template <typename Matrix>
void test(Matrix& matrix, const char* name)
{
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    matrix= 0;
    {
	mtl::mat::inserter<Matrix> ins(matrix);
	for (size_type i= 0; i < matrix.num_rows(); i++)
	    for (size_type j= 0; j < matrix.num_cols(); j++)
		if ((i + j) & 1)
		    ins(i, j) << i + 2*j;
    }

    std::cout << "\n" << name << "\n";
    print_matrix(matrix);

    mtl::mat::transposed_view<Matrix> trans(matrix);
    std::cout << "Transposed" << "\n";
    print_matrix(trans);

    std::cout << "with <<" << "\n"
	      << trans << "\n";

    std::cout << "with << and formatted" << "\n"
	      << with_format(trans, 7, 4) << "\n";

    Matrix square(5, 5);
    square= matrix * trans;

    std::cout << "squared before:\n" << with_format(trans, 4, 2)
	      << "squared in place:\n" << matrix * trans << "\n";

    // Comparison with FP!!!! :-! Make something better eventually
    //MTL_THROW_IF((matrix * trans)[0][1] != 1.0, mtl::runtime_error("Wrong multiplicatin result!"));
}


int main(int, char**)
{
    using namespace mtl;

    dense2D<double>                                dr(5, 7);
    dense2D<double, mat::parameters<col_major> > dc(5, 7);
    morton_dense<double,  morton_mask>             md(5, 7);
    morton_dense<double,  doppled_16_row_mask>     d16r(5, 7);
    compressed2D<double>                           comp(5, 7);

    test(dr, "Dense row major");
    test(dc, "Dense column major");
    test(md, "Morton N-order");
    test(d16r, "Hybrid 16 row-major");
    //test(comp, "compressed2D");

    return 0;
}





