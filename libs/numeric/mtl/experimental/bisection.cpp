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

#include <iostream>
#include <boost/test/minimal.hpp>

#include <boost/numeric/mtl/mtl.hpp>
 
using namespace std;  

template <typename T, typename F>
mtl::maybe<T> inline bisection(T a, T b, const F& f, const T& tau)
{
    if (f(a) * f(b) > 0) 
	return false;
    while (b - a > tau) {
	T m= (a + b) / 2;
	// cout << "f(m) is " << f(m) << '\n';
	if (f(m) * f(a) < 0)
	    b= m;
	else
	    a= m;
    }
    return (a + b) / 2;
}


class caterpillar
{
public:
    caterpillar(double r, double K, double a, double b) : r(r), K(K), a(a), b(b) {}
    
    double operator()(double N) const
    {
	return r * N * (1. - N / K) - a * N * N / (b + N * N);
    }
private:
    double r, K, a, b;
};


template <typename Matrix>
class ev_bisection
{
    typedef typename mtl::Collection<Matrix>::value_type value_type;
public:
    explicit ev_bisection(const Matrix& A) : A(A), D(A), Q(A), R(A)
    {
	assert(num_rows(A) == num_cols(A));
    }

    value_type operator()(const value_type& lambda)
    {
	// D= A - lamdda * identity
	D= -lambda; D+= A;
 
	return lambda; // to make it compile answer is nonsense
    }

private:
    const Matrix &A;
    Matrix D, Q, R;
};


void inline find_eigenvalue(const ev_bisection<mtl::dense2D<double> >& bis, double a, double b)
{
    cout << "Search eigenvalue in interval [" << a << ", " << b << "]\n"; 
    //mtl::maybe<double> lambda= ;
}

int test_main(int argc, char* argv[])
{
#if 0
    caterpillar population1(1.3, 100, 20, 50);
    cout << "Equilibrium of population 1 in interval [0.1, 10] is " << bisection(0.1, 10.0, population1, 0.0001).value() << '\n';
    cout << "Equilibrium of population 1 in interval [10, 20] is " << bisection(10.0, 20.0, population1, 0.0001).value() << '\n';
    cout << "Equilibrium of population 1 in interval [20, 100] is " << bisection(20.0, 100.0, population1, 0.0001).value() << '\n';

    caterpillar population2(2.0, 80, 25, 10);
    cout << "\nEquilibrium of population 2 in interval [0.1, 10] is " << bisection(0.1, 10.0, population2, 0.0001).value() << '\n';
    cout << "Equilibrium of population 2 in interval [10, 20] is " << bisection(10.0, 20.0, population2, 0.0001).value() << '\n';
    cout << "Equilibrium of population 2 in interval [20, 100] is " << bisection(20.0, 100.0, population2, 0.0001).value() << '\n';
#endif    

    mtl::dense2D<double> A(3, 3);
    A= 0.0; A[0][0]= 1.0; A[1][1]= 2.0; A[2][2]= 3;
    
    ev_bisection<mtl::dense2D<double> > bis(A);
    find_eigenvalue(bis, 0, 1.5);

    return 0;
}
