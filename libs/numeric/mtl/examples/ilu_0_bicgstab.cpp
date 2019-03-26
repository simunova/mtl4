#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

using namespace mtl;
using namespace itl;

int main(int, char**)
{
    const int size = 40, N = size * size; 
    typedef compressed2D<double>  matrix_type;

    // Set up a matrix 1,600 x 1,600 with 5-point-stencil
    matrix_type                   A(N, N);
    mat::laplacian_setup(A, size, size);

    // Create an ILU(0) preconditioner
    pc::ilu_0<matrix_type>        P(A);
    
    // Set b such that x == 1 is solution; start with x == 0
    dense_vector<double>          x(N, 1.0), b(N);
    b= A * x; x= 0;
    
    // Termination criterion: r < 1e-6 * b or N iterations
    noisy_iteration<double>       iter(b, 500, 1.e-6);
    
    // Solve Ax == b with left preconditioner P
    bicgstab(A, x, b, P, iter);

    return 0;
}
