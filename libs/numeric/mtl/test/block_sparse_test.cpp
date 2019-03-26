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

// Currently needs -DMTL_DEEP_COPY_CONSTRUCTOR !!!

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace std;
    
    

    // Define a 6x6 sparse matrix in a 3x3 block-sparse
    // Should be 
#if 0 
    typedef mat::parameters<row_major, mtl::index::c_index, fixed::dimensions<2, 2> > parameters1;
    typedef dense2D<double, parameters1>    m_t;
#endif

    // Define a 6x6 sparse matrix in a 3x3 block-sparse
	typedef mtl::dense2D<double>    m_t;
	typedef mtl::compressed2D<m_t>  matrix_t;
    matrix_t                        A(3, 3);
    {
		mtl::mat::inserter<matrix_t> ins(A);

	// First block
	m_t  B(2, 2);
	B= 0.0;
	B[0][0]= 1.0; B[1][1]= 5.0;
	ins(0, 2) << B;

	// Second block
	B=       0.0;
	B[0][1]= 2.0; B[1][0]= 3.0;
	ins(1, 0) << B;

	B= 0.0;
	B[1][0]= 6.0; B[1][1]= 4.0;
	ins(2, 1) << B;
    }
    // cout << "A is " << A << endl; // doesn't works and before it was completely unreadable anyway

    /* Should be something like this:

    [             [ 1 0]] // b1
    [             [ 0 5]] 
    [[ 0 2]             ] // b2
    [[ 3 0]             ]
    [       [ 0 0]      ] // b3
    [       [ 6 4]      ]

    */

    // Access blocks (they are read-only) for sparse matrices
    cout << "The block A(1, 0) is \n" << A(1, 0) << endl;
    cout << "The block A[1][0] is \n" << A[1][0] << endl;

    // Access elements in blocks 
    cout << "In block A(1, 0), the element (0, 1) is " << A(1, 0)(0, 1) << endl;
    cout << "In block A(1, 0), the element [0][1] is " << A(1, 0)[0][1] << endl;
    cout << "In block A[1][0], the element [0][1] is " << A[1][0][0][1] << endl << endl;


	typedef mtl::dense_vector<double> v_t;
	typedef mtl::dense_vector<v_t>    vector_t;
    vector_t                     x(3);

    // x= [[1, 2], [3, 4], [4, 6]]^T
    x[0]= v_t(2, 1.0); x[0][1]= 2.0; // first block of x = [1, 2]^T
    x[1]= v_t(2, 3.0); x[1][1]= 4.0; 
    x[2]= v_t(2, 5.0); x[2][1]= 6.0; 
    cout << "x is " << x << endl; 

    // For y we would only need the vector sizes [[?, ?], [?, ?], [?, ?]]^T
    vector_t                     y(3, v_t(2));

    // Block-sparse matrix * blocked vector !!!
    y= A*x;

    cout << "y after multiplication is " << y << endl
	 << "Should be [[5, 30], [4, 3], [0, 34]]^T." << endl; 

    MTL_THROW_IF(y[1][1] != 3.0, mtl::runtime_error("y[1][1] should be 3.0!\n"));

    return 0;
}
