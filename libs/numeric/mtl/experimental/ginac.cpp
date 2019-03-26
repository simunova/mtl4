#include "boost/numeric/mtl/mtl.hpp"
#include "ginac/ginac.h"
#include <iostream>

using namespace mtl;
using namespace GiNaC;

using mtl::iall;

template <typename Matrix>
void eliminate(Matrix& A, int r, int c)
{
    ex pivot(A[r][c]);
    for (int i= r + 1; i < 4; i++) {
	ex factor(A[i][c] / pivot);
	A[i][iall]-= factor * A[r][iall];
    }
}

template <typename Matrix>
void swap_column(Matrix& A, int c1, int c2)
{
    dense_vector<ex> tmp(clone(A[iall][c1]));
    A[iall][c1]= A[iall][c2];
    A[iall][c2]= tmp;
}

template <typename Matrix, typename Sub>
void substitute(Matrix& A, Sub sub)
{
    for (int i= 0; i < num_rows(A); i++)
	for (int j= 0; j < num_cols(A); j++)
	    A[i][j]= A[i][j].subs(sub);
}

template <typename Matrix>
void trisolve_step(Matrix& A, int r)
{
    for (int i= 3; i > r; i--)
	A[r][4]-= A[r][i] * A[i][4];
    A[r][4]/= A[r][r];
}

//libs:-lcln -lginac
// g++ ginac.cpp -lcln -lginac  -o ginac -I$MTL
int main(int argc, char* argv[])
{
    dense2D< ex > A(4, 5);

    symbol r("r"), s("s");

    A= 1, 0, 2, -3, 2,
	-2, 1, 0, 2, -1,
	-1, 2*r, 6, -5, 4,
	1, 1, 6, r, s;

    std::cout << "A:\n" << A;

    eliminate(A, 0, 0);
    std::cout << "A:\n" << A;

    // eliminate(A, 1, 2);
    // std::cout << "A:\n" << A;

    // swap_column(A, 1, 2);
    // std::cout << "A:\n" << A;

    eliminate(A, 1, 1);
    std::cout << "A:\n" << A;

    eliminate(A, 2, 2);
    std::cout << "A:\n" << A;

    mtl::irange ir(4);
    dense2D<ex> B(clone(A[ir][ir]));

    substitute(A, lst(r == 2, s == 6));
    std::cout << "A:\n" << A;

    // A= 1, 0, 2, 3, 2,
    // 	0, 1, 4, -4, 3,
    // 	0, 0, -8, 0, -4,
    // 	0, 0, 0, 1, 1;

    trisolve_step(A, 3);
    std::cout << "A:\n" << A;

    trisolve_step(A, 2);
    std::cout << "A:\n" << A;   

    trisolve_step(A, 1);
    std::cout << "A:\n" << A;

    trisolve_step(A, 0);
    std::cout << "A:\n" << A;

    substitute(B, lst(r == -7));
    std::cout << "B:\n" << B;

    return 0;
}
