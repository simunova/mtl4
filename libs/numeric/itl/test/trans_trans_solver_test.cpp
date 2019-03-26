#include <iostream>
// #include <boost/mpi.hpp>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


template< typename matrix_type, typename vector_type, typename left_precon_type, 
	  typename right_precon_type, typename iter_type >
void test_solver(const matrix_type& A, vector_type& x, const vector_type& b, 
	   const left_precon_type& L, const right_precon_type& R, iter_type& iter)
{
    itl::bicg(A, x, b, L, iter); // GEHT NICHT!
    itl::bicgstab(A, x, b, L, iter);
    itl::bicgstab_2(A, x, b, L, iter);
    itl::bicgstab_ell(A, x, b, L, R, iter, 3);
    itl::cg(A, x, b, L, R, iter);
    itl::cgs(A, x, b, L, iter);
    itl::gmres(A, x, b, L, R, iter, 30);
    itl::qmr(A, x, b, L, R, iter);
    itl::idr_s(A, x, b, L, R, iter,3);
    itl::tfqmr(A, x, b, L, R, iter);
}

template< typename matrix_type, typename vector_type >
int trans_test(const matrix_type& A, vector_type& x, const vector_type& b)
{
    using namespace mtl;
    using namespace itl;
    
    typedef mat::transposed_view<const matrix_type> trans_matrix_type;
    trans_matrix_type B(A);
    
    // Create an ILU(0) preconditioner
    pc::ilu_0<matrix_type>        	P0(A);
    pc::ilu_0<trans_matrix_type>	P1(B);
    
    pc::identity<matrix_type>		Id0(A);
    pc::identity<trans_matrix_type> 	Id1(B);
  
    // Termination criterion: r < 1e-6 * b or N iterations
    noisy_iteration<double>       iter(b, 500, 1.e-6);
    
    test_solver(A, x, b, P0, Id0, iter);
    test_solver(A, x, b, P1, Id1, iter); // GEHT NICHT!
    test_solver(B, x, b, P0, Id0, iter);
    test_solver(B, x, b, P1, Id1, iter); // GEHT NICHT!

    return 0;
}

int main(// int argc, char* argv[]
	 ) 
{
    using namespace mtl;

    // mtl::par::environment env(argc, argv);

    const int size = 10, N = size * size;
    // typedef mat::distributed<compressed2D<double> >  matrix_type;
    typedef compressed2D<double>   matrix_type;
    matrix_type                                         A;
    laplacian_setup(A, size, size);
    
    mtl::dense_vector<double>                           x(N, 1.0), b;
    
    b= A * x;
    x= 0;
    
    int error_code = trans_test(A, x, b);
    
    return error_code;
}
