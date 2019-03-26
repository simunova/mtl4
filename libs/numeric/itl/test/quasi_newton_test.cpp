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
// #include <boost/test/minimal.hpp>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

using namespace std;  
   
template <typename Vector>
struct f_test
{
    f_test() : s(3)
    {
	mtl::vec::inserter<Vector> ins(s);
	ins[0] << 1; ins[1] << 2; ins[2] << 2; 
    }

    typename mtl::Collection<Vector>::value_type 
    operator() (const Vector& x) const
    {
	return dot(s, Vector(ele_prod(x, x)));
    }
    Vector s;
};


template <typename Vector>
struct grad_f_test
{
    grad_f_test() : s(3)
    {
	mtl::vec::inserter<Vector> ins(s);
	ins[0] << 2; ins[1] << 4; ins[2] << 4; 
    }

    Vector operator() (const Vector& x) const {  return Vector(ele_prod(s, x)); }
    Vector s;
};


int main(int, char**)
{
    using namespace mtl;
    using mtl::io::tout;

    mtl::dense_vector<double>       x(3, 8);
    tout << "x= " << x << "\n";

    grad_f_test<mtl::dense_vector<double> > grad_f;
    f_test<mtl::dense_vector<double> >      f;
        
    itl::cyclic_iteration<double>  iter(grad_f(x), 1000, 0, 1e-4, 100);
    itl::cyclic_iteration<double> iter1(grad_f(x), 1000, 0, 1e-4, 100);
    itl::cyclic_iteration<double> iter2(grad_f(x), 1000, 0, 1e-4, 100);
    itl::cyclic_iteration<double> iter3(grad_f(x), 1000, 0, 1e-4, 100);
    itl::cyclic_iteration<double> iter4(grad_f(x), 1000, 0, 1e-4, 100);
    
    quasi_newton(x, f, grad_f, itl::wolf<>(), itl::bfgs(), iter);
    iter.error_code();    

   // tout << "x= " << x << "\n";
    tout << "grad_f(x)= " << grad_f(x) << "\n\n";
    if (two_norm(x) > 10 * iter.atol()) throw "x should be 0.";
    x= 8;
    quasi_newton(x, f, grad_f, itl::wolf<>(), itl::dfp(), iter1);
    iter1.error_code();    

   // tout << "dfp x= " << x << "\n";
    tout << "grad_f(x)= " << grad_f(x) << "\n\n";
    if (two_norm(x) > 10 * iter1.atol()) throw "x should be 0.";
    
    x= 8;
    quasi_newton(x, f, grad_f, itl::wolf<>(), itl::broyden(), iter2);
    iter2.error_code();    

    tout << "broyden x= " << x << "\n";
    tout << "grad_f(x)= " << grad_f(x) << "\n\n";
    if (two_norm(x) > 10 * iter2.atol()) throw "x should be 0.";
    
 #if 0  //bad condition on some compiler
    x= 8;
    quasi_newton(x, f, grad_f, itl::wolf<>(), itl::sr1(), iter3);
    iter3.error_code();    

    tout << "sr1 x= " << x << "\n";
    tout << "grad_f(x)= " << grad_f(x) << "\n";
    if (two_norm(x) > 10 * iter3.atol())
	throw "x should be 0.";
#endif     

    x= 8;
    quasi_newton(x, f, grad_f, itl::wolf<>(), itl::psb(), iter4);
    iter4.error_code();

    tout << "psb x= " << x << "\n";
    tout << "grad_f(x)= " << grad_f(x) << "\n\n";
    if (two_norm(x) > 10 * iter4.atol()) throw "x should be 0.";


    return 0;
}
 














