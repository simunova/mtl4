#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl;
    
    double array[2][3] = {{8, 9, 11},
			  {2, 6, 17}};
    dense2D<double>      A(array);

    trans(A)[1][0]= 7;    // Modify A[0][1]

    std::cout << "A is\n" << A << "\n";

    const dense2D<double> B(array);
 
    //  trans(B)[1][0]= 7;   // Error: transposed of B is also constant
    std::cout << "trans(B) is\n" << trans(B) << "\n";  // Read access is of course possible

    return 0;
}
