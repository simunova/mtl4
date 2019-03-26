// File: multi_vector.cpp

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;

    typedef dense_vector<double>   Vector;

    Vector                      v(2, 3.4), w(3, 2.5);
    mtl::multi_vector<Vector> 	A(2, 3);    
    dense2D<double>		B(2,2), C(3,2), D(3,3);

    // Initialize matrices
    A= 3.0; B= 4.0; C= 5.0; D= 6.0;

    // vector= multi_vector * vector
    v= A * w;

    // vector= transposed multi_vector * vector
    w= trans(A) * v;

    // vector= matrix * vector 
    v= B * A.vector(1);		

    // vector= matrix * vector
    A.vector(0)= B * A.vector(1);	

    // Orthogonalize multi_vector
    orth(A);

    return 0;
}

