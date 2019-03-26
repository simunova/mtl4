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

// #define MTL_VERBOSE_TEST

#include <iostream>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

const int size = 3, N = size * size; 

int max_iter= 0;

template <typename Matrix, typename Solver>
void test1(const char* solver_name, const char* matrix_name)
{
    mtl::io::tout << solver_name << " on " << matrix_name << std::endl;;

    Matrix                             A;
    laplacian_setup(A, size, size);
       
    mtl::dense_vector<double>          x(N, 1.0), b(N);    
    b = A * x;
    x= 0;
    
    Solver s(A);
    s.iteration_ref().set_max_iterations(N+1);

#ifdef MTL_VERBOSE_TEST
    s.iteration_ref().set_quite(false);
    s.iteration_ref().set_cycle(1);
#else
    s.iteration_ref().suppress_resume(true);
#endif

    s(x, b);
    int i= s.iteration().iterations();
    if (i > max_iter)
	max_iter= i;
    // MTL_THROW_IF(size == 3 && i > 9, mtl::runtime_error("Too many iterations in solver"));
    if (size == 3 && i > 9)
	std::cout << "Solver \"" << solver_name << "\" converges slowly!\n";
}

// same without iter
template <typename Matrix, typename Solver>
void test1a(const char* solver_name, const char* matrix_name)
{
    mtl::io::tout << solver_name << " on " << matrix_name << std::endl;;

    Matrix                             A;
    laplacian_setup(A, size, size);
       
    mtl::dense_vector<double>          x(N, 1.0), b(N);    
    b = A * x;
    x= 0;
    
    Solver s(A);
    s.step(x, b);
}



template <typename Matrix>
int test2(const char* matrix_name)
{
    test1<Matrix, itl::cg_solver<Matrix> >("CG", matrix_name);
    test1<Matrix, itl::cg_solver<Matrix, itl::pc::ilu_0<Matrix> > >("CG with ILU_0", matrix_name);
    test1<Matrix, itl::cgs_solver<Matrix> >("CGS", matrix_name);
    test1<Matrix, itl::cgs_solver<Matrix, itl::pc::ilu_0<Matrix> > >("CGS with ILU_0", matrix_name);

    test1<Matrix, itl::bicg_solver<Matrix> >("BiCG", matrix_name);
    test1<Matrix, itl::bicg_solver<Matrix, itl::pc::ilu_0<Matrix> > >("BiCG with ILU_0", matrix_name);
    test1<Matrix, itl::bicgstab_solver<Matrix> >("BiCGStab", matrix_name);
    test1<Matrix, itl::bicgstab_solver<Matrix, itl::pc::ilu_0<Matrix> > >("BiCGStab with ILU_0", matrix_name);
    test1<Matrix, itl::bicgstab_2_solver<Matrix> >("BiCGStab(2)", matrix_name);
    test1<Matrix, itl::bicgstab_2_solver<Matrix, itl::pc::ilu_0<Matrix> > >("BiCGStab(2) with ILU_0", matrix_name);
    test1<Matrix, itl::bicgstab_ell_solver<Matrix> >("BiCGStab(ell)", matrix_name);
    test1<Matrix, itl::bicgstab_ell_solver<Matrix, itl::pc::ilu_0<Matrix> > >("BiCGStab(ell) with ILU_0", matrix_name);
    test1<Matrix, itl::gmres_solver<Matrix> >("GMRES", matrix_name);
    test1<Matrix, itl::gmres_solver<Matrix, itl::pc::ilu_0<Matrix> > >("GMRES with ILU_0", matrix_name);
    test1<Matrix, itl::gmres_solver<Matrix, itl::pc::ilu_0<Matrix>, itl::pc::ilu_0<Matrix> > >("GMRES with ILU_0 from left and right", matrix_name);

    test1<Matrix, itl::qmr_solver<Matrix> >("QMR", matrix_name);
    test1<Matrix, itl::qmr_solver<Matrix, itl::pc::ilu_0<Matrix> > >("QMR with ILU_0", matrix_name);
    test1<Matrix, itl::tfqmr_solver<Matrix> >("TFQMR", matrix_name);
    test1<Matrix, itl::tfqmr_solver<Matrix, itl::pc::ilu_0<Matrix> > >("TFQMR with ILU_0", matrix_name);

    test1<Matrix, itl::idr_s_solver<Matrix> >("IDR(s)", matrix_name);
    test1<Matrix, itl::idr_s_solver<Matrix, itl::pc::ilu_0<Matrix> > >("IDR(s) with ILU_0", matrix_name);

    typedef itl::cg_solver<Matrix>                           cg1_type;
    typedef itl::cg_solver<Matrix, itl::pc::ilu_0<Matrix> >  cg2_type;
    
    test1a<Matrix, itl::repeating_solver<cg1_type, 3, true> >("3 iterations CG", matrix_name);
    test1a<Matrix, itl::repeating_solver<cg2_type, 3, true> >("3 iterations CG with ILU_0", matrix_name);

    return 0;
}

int main()
{
    test2<mtl::compressed2D<double> >("Compressed2D");
    //test<mtl::sparse_banded<double> >();
    // std::cout << "Largest iteration number is " << max_iter << '\n';

    return 0;
}
