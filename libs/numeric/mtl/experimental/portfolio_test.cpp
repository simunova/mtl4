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
#include <utility>
#include <cmath>
#include <boost/test/minimal.hpp>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

using namespace std;  
using namespace mtl;  
   

template <typename Vector, typename Matrix>
class f_ftor
{
    typedef typename Collection<Vector>::value_type value_type;

  public:
    /// Arguments: each stock's ROI, the aimed ROI, the Langrange factor and the covariance
    f_ftor(const Vector& rv, value_type r, value_type lg, const Matrix& S) 
      : rv(rv), r(r), lg(lg), S(S) {}

    value_type operator()(const Vector& pi) const
    {
	Vector S_pi(S * pi);
	return lg * sq(sum(pi) - 1) + lg * sq(dot(pi, rv) - r) + trans(pi) * S_pi;
    }
    value_type sq(value_type x) const { return x * x; }

  private:
    Vector      rv;
    value_type  r, lg;
    Matrix      S;
};


template <typename Vector, typename Matrix>
class grad_f_ftor
{
    typedef typename Collection<Vector>::value_type value_type;

  public:
    /// Arguments: each stock's ROI, the aimed ROI, the Langrange factor and the covariance
    grad_f_ftor(const Vector& rv, value_type r, value_type lg, const Matrix& S) 
      : rv(rv), onev(size(rv), 1), r(r), lg(lg), S(S) {}

    Vector operator()(const Vector& pi) const
    {
	value_type f1= 2.0 * lg * (sum(pi) - 1), f2= 2.0 * lg * (dot(pi, rv) - r);
	return Vector(f1 * onev + f2 * rv + 2.0 * S * pi);
    }

  private:
    Vector      rv, onev;
    value_type  r, lg;
    Matrix      S;
};


template <typename Vector, typename Matrix>
class portfolio_optimizer
{
    typedef typename Collection<Vector>::value_type value_type;
    
  public:
    /// Arguments: each stock's ROI, the aimed ROI and the covariance
    portfolio_optimizer(const Vector& rv, value_type r, const Matrix& S) 
      : s(size(rv) + 2), A(s, s), b(s, value_type(0))
    {
	unsigned s2= size(rv);
	A[irange(s2)][irange(s2)]= S; A[irange(s2)][s2]= Vector(s2, 1); A[irange(s2)][s2+1]= rv;
	A[s2][irange(s2)]= trans(Vector(s2, 1)); A[irange(s2, s)][irange(s2, s)]= 0;
	A[s2+1][irange(s2)]= trans(rv);

	b[s2]= 1; b[s2+1]= r;
	cout << "A is\n" << A << "\nb is " << b << '\n';
    }

    Vector operator()() const {	return clone(lu_solve(A, b)[irange(s-2)]); }

  private:

    unsigned    s;
    Matrix      A;
    Vector      b;
 };



int test_main(int, char**)
{
    using namespace mtl;

    dense_vector<double>       pi(4, 0.25), rv(4);
    rv= 1.03, 1.14, 1.05, 1.08;

    dense2D<double>            S(4, 4);
#if 1
    S= 1.0, 0.1, 0.3, -0.2,
	0.1, 1., -0.4, 0.7,
	0.3, -0.4, 1., 0.4,
	-0.2, 0.7, 0.4, 1.;
#else
    S= 1.0, 0.1, 0.3, 0.2,
	0.1, 1., 0.4, 0.7,
	0.3, 0.4, 1., 0.4,
	0.2, 0.7, 0.4, 1.;
#endif

    const double lagrange= 10000.0;
    f_ftor< dense_vector<double>, dense2D<double> >  f(rv, 1.09, lagrange, S);
    cout << "f(pi) is " << f(pi) << '\n';

    grad_f_ftor< dense_vector<double>, dense2D<double> >  grad_f(rv, 1.09, lagrange, S);
    cout << "grad_f(pi) is " << grad_f(pi) << '\n';
    
    portfolio_optimizer< dense_vector<double>, dense2D<double> >  opt(rv, 1.09, S);
    dense_vector<double>       pi_opt(opt()); 
    cout << "Optimal portfolio is " << pi_opt << "\n";

    std::cout<< "Sum of pi is " << sum(pi_opt) << "\n";
    std::cout<< "Overall ROI is " << dot(pi_opt, rv) << "\n";
    std::cout<< "Variance is " << dot(pi_opt, dense_vector<double>(S * pi_opt)) << "\n";

#if 1
    std::cout<< "Sum of pi is " << sum(pi) << "\n";
    std::cout<< "Overall ROI is " << dot(pi, rv) << "\n";
    std::cout<< "Variance is " << dot(pi, dense_vector<double>(S * pi)) << "\n";

    itl::cyclic_iteration<double> iter(grad_f(pi), 1000, 0, 1e-5, 10);
    quasi_newton(pi, f, grad_f, itl::armijo<>(), itl::sr1(), iter);
    iter.error_code();    

    std::cout<< "pi= " << pi << "\n";
    std::cout<< "grad_f(pi)= " << grad_f(pi) << "\n";

    std::cout<< "Sum of pi is " << sum(pi) << "\n";
    std::cout<< "Overall ROI is " << dot(pi, rv) << "\n";
    std::cout<< "Variance is " << dot(pi, dense_vector<double>(S * pi)) << "\n";
#endif

    return 0;
}
 














