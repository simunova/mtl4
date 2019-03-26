// Filename: matrix_free_pcg.cpp (part of MTL4)

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

struct poisson2D_diagonal_pc
{
    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& x, VectorOut& y) const
    {
	for (std::size_t i= 0; i < size(x); i++)
	    y[i]= x[i] * 0.25;
    }

    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& x, VectorOut& y) const
    {
	for (std::size_t i= 0; i < size(x); i++)
	    y[i]= x[i] * 0.25;
    }
};

template <typename Vector>
itl::pc::solver<poisson2D_diagonal_pc, Vector, false>
inline solve(const poisson2D_diagonal_pc& P, const Vector& x)
{
    return itl::pc::solver<poisson2D_diagonal_pc, Vector, false>(P, x);
}

template <typename Vector>
itl::pc::solver<poisson2D_diagonal_pc, Vector, true>
inline adjoint_solve(const poisson2D_diagonal_pc& P, const Vector& x)
{
    return itl::pc::solver<poisson2D_diagonal_pc, Vector, true>(P, x);
}


int main(int, char**)
{
    // For a more realistic example set size to 1000 or larger
    const int size = 10, N = size * size;

    typedef mtl::mat::poisson2D_dirichlet  matrix_type;
    matrix_type                               A(size, size);
    poisson2D_diagonal_pc                     P;

    mtl::dense_vector<double>                 x(N, 1.0), b(N);

    b = A * x;
    x= 0;
    itl::cyclic_iteration<double>             iter(b, 100, 1.e-11, 0.0, 5);
    cg(A, x, b, P, P, iter);

    return 0;
}

