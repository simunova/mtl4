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

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/timer.hpp>

template <typename Matrix>
void setup(Matrix& A)
{
    const int n= num_rows(A);
    A= 1.0;
    mtl::mat::inserter<Matrix, mtl::update_plus<double> > ins(A);

    for (int i= 0; i < 30 * n; i++) {
	int r= rand()%n, c= rand()%n;
	ins[r][c] << -1;
	ins[r][r] << 1;
    }
}



int main(int argc, char* argv[])
{
    // For a more realistic example set sz to 1000 or larger
    int size = 3;
    if (argc > 1) size = atoi(argv[1]);
    int N = size * size; 

    typedef mtl::compressed2D<double>  matrix_type;
    typedef mtl::dense_vector<double>  vector_type;
    mtl::compressed2D<double>          A(N, N);
    setup(A);

    boost::timer                       fac_timer;
    itl::pc::ilu_0<matrix_type>        P(A);
    std::cout << "Factorization took " << fac_timer.elapsed() << "s\n";
    
    mtl::dense_vector<double> x(N, 3.0), x2(N), Px(N), x3(N), x4(N), x5(N);

    matrix_type L(P.get_L()), U(P.get_U());

    x2= strict_upper(U) * x;
    for (int i= 0; i < N; i++)
	x2[i]+= 1. / U[i][i] * x[i];

    Px= L * x2 + x2;

    x4= unit_lower_trisolve(L, Px);
    if (two_norm(vector_type(x4 - x2)) > 0.01) throw "Error in unit_lower_trisolve.";

    x5= inverse_upper_trisolve(U, x4);
    if (two_norm(vector_type(x5 - x)) > 0.01) throw "Error in inverse_upper_trisolve.";


    boost::timer                       solve_timer;
    x3= solve(P, Px);
    std::cout << "Solving took " << solve_timer.elapsed() << "s\n";

    if (two_norm(vector_type(x3 - x)) > 0.01) throw "Error in solve.";

    return 0;
}
