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

void test(int n, int m, int order)
{
    assert(n%3 == 0);
    mtl::dense2D<double>           block(3, m);
    mtl::dense_vector<std::size_t> rows(3), columns(m);

    for (int i= 0; i < 3; i++) rows[i]= i;

    switch(order) {
    case 0: for (int i= 0; i < m; i++) columns[i]= i; break;
    case 1: for (int i= 0; i < m; i++) columns[i]= m-i-1; break;
    case 2: for (int i= 0; i < m; i++) columns[i]= (m/2+i)%m; break;
    }
    columns[3]= columns[1]; columns[2]= columns[0];

    for (int i= 0; i < m; i++) {
	block[0][i]= 10+i;
	block[1][i]= 20+i;
	block[2][i]= 30+i;
    }

    mtl::compressed2D<double> A(n, n);
    {
        mtl::mat::inserter<mtl::compressed2D<double>, mtl::update_plus<double> > ins(A, 5);
	for (unsigned r= 0; r < num_rows(A); r+= 3) {
	    ins << element_matrix(block, rows, columns);
	    ins << element_matrix(block, rows, columns);
	    for (int i= 0; i < 3; i++) rows[i]+= 3;
	    for (int i= 0; i < m; i++) ++columns[i];
	}
    }
    if (n < 11) cout << "A is \n" << with_format(A, 4, 3);
    switch(order) {
    case 0:
	MTL_THROW_IF(A[1][2] != 0, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF(A[4][2] != 88, mtl::runtime_error("Wrong value"));break;
    case 1:
	MTL_THROW_IF(A[1][2] != 0, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF(A[4][2] != 0, mtl::runtime_error("Wrong value"));	break;
    case 2:
	MTL_THROW_IF(A[1][2] != 84, mtl::runtime_error("Wrong value"));
	MTL_THROW_IF(A[4][2] != 48, mtl::runtime_error("Wrong value"));	break;
    }

    double array[4]= {1.0, -.4, -0.5, 2.0};
    int v0[]= {6, 7}, v1[]= {7, 8};
    {
        mtl::mat::inserter<mtl::compressed2D<double>, mtl::update_plus<double> > inserter(A, 5);
	inserter << mtl::element_array(mtl::dense2D<double>(mtl::size(v0), mtl::size(v1), array), v0, v1);
    }
    if (n < 11) cout << "A is \n" << with_format(A, 4, 3);
    MTL_THROW_IF(A[7][7] != -0.5, mtl::runtime_error("Wrong value"));

}

 
int main(int, char**)
{
    test(9, 5, 0);
    test(9, 5, 1);
    test(9, 5, 2);

    return 0;
}
