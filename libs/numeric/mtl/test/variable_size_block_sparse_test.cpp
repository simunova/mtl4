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


template <typename T>
struct as {};

template <>
struct as<int> 
{
	typedef mtl::ashape::scal  type;
};

template <typename V>
struct as<mtl::compressed2D<V> >
{
	typedef mtl::ashape::mat<typename as<V>::type>   type;
};

template <typename V>
struct as<mtl::dense2D<V> >
{
	typedef mtl::ashape::mat<typename as<V>::type>   type;
};

int main(int, char**)
{
    using namespace std;

#if 0  
    cout << typeid(as<int>::type).name() << endl;
    cout << typeid(as<mtl::compressed2D<int> >::type).name() << endl;
    cout << typeid(as<mtl::compressed2D<mtl::dense2D<int> > >::type).name() << endl;

    cout << typeid(mtl::ashape::ashape<int>::type).name() << endl;
    cout << typeid(mtl::ashape::ashape<mtl::compressed2D<int> >::type).name() << endl;
    cout << typeid(mtl::ashape::ashape<typename mtl::compressed2D<mtl::dense2D<int> >::value_type >::type).name() << endl;
#endif

    // Define a 6x5 sparse matrix in a 3x3 block-sparse
    typedef mtl::dense2D<double>    m_t;
    typedef mtl::compressed2D<m_t>  matrix_t;
    matrix_t                        A(3, 3);
    {
	mtl::mat::inserter<matrix_t> ins(A);

	// First block
	m_t  b1(1, 1);
	b1(0, 0)= 1.0;
	ins(0, 2) << b1;

	// Second block
	m_t  b2(2, 3);
	b2=       0.0;
	b2[0][1]= 2.0;
	b2[1][2]= 3.0;
	ins(1, 0) << b2;

	m_t b3(3, 1);
	b3= 0.0;
	b3[1][0]= 4.0;
	ins(2, 1) << b3;
    }
    // cout << "A is " << A << endl; // doesn't works and before it was completely unreadable anyway

    /* Should be something like this:

       [                   [ 1]] // b1
       [[  0   2   0]          ] // b2
       [[  0   0   3]          ]
       [             [  0]     ] // b3
       [             [  4]     ]
       [             [  0]     ]

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
    vector_t                          x(3), y(3);

    // x= [[0, 5, 3], [1], [8]]^T
    x[0]= v_t(3, 0.0); x[0][1]= 5.0; x[0][2]= 3.0; // first block of x = [0, 5, 3]^T
    x[1]= v_t(1, 1.0);                             // second block of x 
    x[2]= v_t(1, 8.0);                             // third block

    cout << "x is " << x << endl; 

    // For y we would only need the vector sizes [[?], [?, ?], [?, ?, ?]]^T
    // To avoid valgrind complains we set to 0
    y[0]= v_t(1, 0.0); y[1]= v_t(2, 0.0); y[2]= v_t(3, 0.0);

    cout << "y is " << y << endl; 

    // Block-sparse matrix * blocked vector !!!
    y= A*x;

    cout << "y after multiplication is " << y << endl
	 << "Should be [[8], [10, 9], [0, 4, 0]]^T." << endl; 

    return 0;
}
