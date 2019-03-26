#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main()
{
  const int size = 3, N = size * size;
  mtl::compressed2D<double>          A(N, N);
 
  laplacian_setup(A, size, size);
  std::cout<< "Laplacian_setup=\n" << A << "\n";
  
  diagonal_setup(A, 2.0);
  std::cout<< "diagonal_setup=\n" << A << "\n";
  
  mtl::dense2D<double>		     B(N, N);   //to expensive for sparse Matrix
  hessian_setup(B, 4.0);
  std::cout<< "hessian_setup=\n" << B << "\n";
    
  return 0;
}
