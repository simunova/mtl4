// Filename: umfpack_solve_example.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;  

int main(int, char**)
{
#ifdef MTL_HAS_UMFPACK
    typedef mtl::compressed2D<double> matrix_type;

    matrix_type A(5, 5);
    A= 2.,  3.,  0.,  0.,  0.,
       3.,  0.,  4.,  0.,  6.,
       0., -1., -3.,  2.,  0.,
       0.,  0.,  1.,  0.,  0.,
       0.,  4.,  2.,  0.,  1.;
    crop(A);

    mtl::dense_vector<double>   x(5), b(5);
    b= 8., 45., -3., 3., 19.;
    mtl::dense_vector<double>   b2(2 * b);
    cout << "A = \n" << A << "b = " << b << "\n";

    // Factorize and solve
    umfpack_solve(A, x, b);
    cout << "\nA \\ b using umfpack_solve = " << x << "\n";
    
    // Define a solver object by internally factorizing A
    mtl::mat::umfpack::solver<matrix_type> solver(A);

    // Solve A * x == b and b2 with the solver object
    solver(x, b);
    solver(x, b2);

    // Change one or more matrix entries while keeping the sparsity pattern
    A.lvalue(1, 2)= 5.0;

    // Compute a new factorization (relying on unchanged sparsity)
    solver.update_numeric();
    
    // If we change b accordingly we will get the same result
    b[1]= 48;
    solver(x, b);
    cout << "\nA \\ b after numeric update = " << x << "\n";

    // Change matrix's values and sparsity
    {
	mtl::mat::inserter<matrix_type> ins(A);
	ins[3][4] << 2.;	
    }
    cout << "\nA is now = \n" << A << "\n";

    // Perform a completely new factorization
    solver.update();

    b[3]= 13.;
    int status= solver(x, b);
    cout << "A \\ b after (complete) update = " << x << ", status is " << status << "\n";

#endif
    return 0;
}
