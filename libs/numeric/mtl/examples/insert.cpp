// File: insert.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace mtl;

template <typename Matrix>
void fill(Matrix& m)
{
    // Matrices are not initialized by default
    m= 0.0;

    // Create inserter for matrix m
    mat::inserter<Matrix> ins(m);
    
    // Insert value in m[0][0]
    ins[0][0] << 2.0;
    ins[1][2] << 0.5;
    ins[2][1] << 3.0;

    // Destructor of ins sets final state of m
}

template <typename Matrix>
void modify(Matrix& m)
{
    // Type of m's elements
    typedef typename Collection<Matrix>::value_type value_type;

    // Create inserter for matrix m
    // Existing values are not overwritten but inserted
    mat::inserter<Matrix, update_plus<value_type> > ins(m, 3);
    
    // Increment value in m[0][0]
    ins[0][0] << 1.0;

    // Elements that doesn't exist (in sparse matrices) are inserted
    ins[1][1] << 2.5;
    ins[2][1] << 1.0;
    ins[2][2] << 4.0;

    // Destructor of ins sets final state of m
}

int main(int, char**)
{
    // Matrices of different types
    compressed2D<double>              A(3, 3);
    dense2D<double>                   B(3, 3);
    morton_dense<float, morton_mask>  C(3, 3);

    // Fill the matrices generically
    fill(A); fill(B); fill(C);
    std::cout << "A is \n" << A << "\nB is \n" << B << "\nC is \n" << C;

    // Modify the matrices generically
    modify(A); modify(B); modify(C);
    std::cout << "\n\nAfter modification:\nA is \n" << A 
	      << "\nB is \n" << B << "\nC is \n" << C;

    return 0;
}

