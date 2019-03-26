// File: element_matrix.cpp

#include <iostream>
#include <vector>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

template <typename Matrix>
void fill(Matrix& m)
{
    // Matrices are not initialized by default
    m= 0.0;

    // Type of m's elements
    typedef typename Collection<Matrix>::value_type value_type;

    // Create inserter for matrix m
    // Existing values are not overwritten but inserted
    mat::inserter<Matrix, update_plus<value_type> > ins(m, 3);
    
    // Define element matrix (array)
    double m1[2][2]= {{1.0, -.4}, {-0.5, 2.0}}; 

    // Corresponding indices of the elements
    std::vector<int> v1(2);
    v1[0]= 1; v1[1]= 3;

    // Insert element matrix
    ins << element_array(m1, v1);

    // Insert same array with different indices
    v1[0]= 0; v1[1]= 2;
    ins << element_array(m1, v1);

    // Use element matrix type with dynamic size
    dense2D<double> m2(2, 3);
    m2[0][0]= 1; m2[0][1]= 0.2; m2[0][2]= 0.1; 
    m2[1][0]= 2; m2[1][1]= 1.2; m2[1][2]= 1.1;

    // Vector for column indices 
    dense_vector<int> v2(3);
    // Indices can be out of order
    v2[0]= 4; v2[1]= 1; v2[2]= 3;

    // Use element_matrix and separate vectors for row and column indices
    ins << element_matrix(m2, v1, v2);
}

int main(int, char**)
{
    // Matrices of different types
    compressed2D<double>              A(5, 5);
    dense2D<double>                   B(5, 5);
    morton_dense<float, morton_mask>  C(5, 5);

    // Fill the matrices generically
    fill(A); fill(B); fill(C);
    std::cout << "A is \n" << with_format(A, 4, 3) 
	      << "\nB is \n" << with_format(B, 4, 3)
	      << "\nC is \n" << with_format(C, 4, 3);

    return 0;
}

