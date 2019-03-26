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

#undef MTL_ASSERT_FOR_THROW // Don't wont to abort program when one solver fail

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// #include <boost/test/minimal.hpp>

#define MTL_RUN_SOLVER( name, solver, args)				\
    {									\
	std::cout << "\n\n" << name << "\n";				\
	xs= 0.01;							\
        int codes= 0;							\
        itl::cyclic_iteration<double> iters(bs, 150, 1.e-4, 0.0, 10);     \
        try {								\
	    codes= solver args;						\
	} catch (const itl::unexpected_orthogonality&) {		\
	    std::cerr << "Unexpected_Orthogonality in solver!\n";	\
	}								\
        if (codes != 0) {						\
	    std::cerr << "Solver doesn't converge!!!\n";		\
	    failed++;							\
        } else								\
	    succeed++;							\
    }
    
 
int main(int, char**)
{
    // For a more realistic example set size to 1000 or larger
    const int size = 4, N = size * size;
    int       succeed= 0, failed= 0;
    
    typedef mtl::compressed2D<double>  matrix_s_type;
    matrix_s_type                                           As(N, N);
    laplacian_setup(As, size, size);
    mtl::dense_vector<double>                               xs(N, 1.0), bs(N);
    bs= As * xs;
    itl::pc::identity<matrix_s_type>                        I(As);
    itl::pc::ic_0<matrix_s_type>                            IC(As);
    itl::pc::ilu_0<matrix_s_type>                           ILU(As);
    const unsigned                                          ell= 6, restart= 8, s= ell;

    MTL_RUN_SOLVER("Bi-Conjugate Gradient", bicg, (As, xs, bs, I, iters));
    MTL_RUN_SOLVER("Bi-Conjugate Gradient Stabilized", bicgstab, (As, xs, bs, ILU, iters));
    MTL_RUN_SOLVER("Bi-Conjugate Gradient Stabilized(2)", bicgstab_2, (As, xs, bs, ILU, iters));
    MTL_RUN_SOLVER("Bi-Conjugate Gradient Stabilized(ell)", bicgstab_ell, (As, xs, bs, ILU, I, iters, ell));
    MTL_RUN_SOLVER("Conjugate Gradient", cg, (As, xs, bs, IC, iters));
    MTL_RUN_SOLVER("Conjugate Gradient Squared", cgs, (As, xs, bs, ILU, iters));
//     MTL_RUN_SOLVER("Generalized Minimal Residual method (without restart)", gmres_full, (As, xs, bs, I, I, iters, size));  //only N iterations
    MTL_RUN_SOLVER("Generalized Minimal Residual method with restart", gmres, (As, xs, bs, I, I, iters, restart));
    MTL_RUN_SOLVER("Induced Dimension Reduction on s dimensions (IDR(s))", idr_s, (As, xs, bs, ILU, I, iters, s));
    MTL_RUN_SOLVER("Quasi-minimal residual", qmr, (As, xs, bs, ILU, I, iters));
    MTL_RUN_SOLVER("Transposed-free Quasi-minimal residual", tfqmr, (As, xs, bs, ILU, I, iters));
    std::cout << succeed << " solvers succeeded and " << failed << " solvers failed.\n";
  
    return 0;
}
