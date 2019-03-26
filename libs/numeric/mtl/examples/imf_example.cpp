// Filename: imf_example.cpp (part of MTL4)

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

int main(int, char** argv)
{
    typedef double value_type;
  
    std::string program_dir= mtl::io::directory_name(argv[0]),
  	        matrix_file= mtl::io::join(program_dir, "../../mtl/test/matrix_market/square3.mtx");
    //define and read element structure
    mtl::mat::element_structure<value_type> A;
    read_el_matrix(matrix_file, A);
    
    int size= int( num_cols(A) );
    mtl::dense_vector<value_type>    x(size, 1), b( A * x );
    
    //assemble sparse matrix from element structure (unnecessary)
    mtl::compressed2D<value_type> B;
    assemble_compressed(A, B);
    
    //create imf preconditioner with 3 levels of fill-in
    itl::pc::imf_preconditioner<value_type> precond(A, 3);
    
    itl::cyclic_iteration<value_type>    iterA(b, size, 1.e-8, 0.0, 5),
					  iterB(b, size, 1.e-8, 0.0, 5);
    x= 0;
    bicgstab(A, x, b, precond, iterA);  //solve A*x=b  with element_structure A
    x= 0;
    bicgstab(B, x, b, precond, iterB);  //solve A*x=b  with assembled sparse matrix B
    
    return 0;
}
