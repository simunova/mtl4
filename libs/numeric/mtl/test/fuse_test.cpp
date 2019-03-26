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
#include <boost/numeric/mtl/mtl.hpp>

inline bool close(double x, double y) 
{ 
    using std::abs;
    return abs(x - y) < 0.001;
}

int main(int , char**)
{
    using namespace std;
    using mtl::lazy;
    
    double                d, rho, alpha= 7.8, beta, gamma;
    std::complex<double>  z;

    mtl::dense_vector<double> v(6, 1.0), w(6), r(6, 6.0), q(6, 2.0), x(6);
    mtl::dense2D<double>      A(6, 6);
    A= 2.0;
    mtl::compressed2D<double>      B(6, 6);
    B= 2.0;

    (lazy(w)= A * v) || (lazy(d) = lazy_dot(w, v));
    cout << "w = " << w << ", d = " << d << "\n";
    MTL_THROW_IF(!close(d, 12), mtl::runtime_error("wrong dot"));    

    (lazy(w)= B * v) || (lazy(d) = lazy_dot(w, v));
    cout << "w = " << w << ", d = " << d << "\n";
    MTL_THROW_IF(!close(d, 12), mtl::runtime_error("wrong dot"));

    (lazy(r)-= alpha * q) || (lazy(rho)= lazy_unary_dot(r)); 
    cout << "r = " << r << ", rho = " << rho << "\n";
    MTL_THROW_IF(!close(rho, 552.96), mtl::runtime_error("wrong unary_dot"));

    (lazy(x)= 7.0) || (lazy(beta)= lazy_unary_dot(x)); 
    cout << "x = " << x << ", beta = " << beta << "\n";
    MTL_THROW_IF(!close(beta, 294), mtl::runtime_error("wrong unary_dot"));
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_one_norm(x)); 
    cout << "x = " << x << ", beta = " << beta << "\n";
    MTL_THROW_IF(!close(beta, 42), mtl::runtime_error("wrong one_norm"));
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_two_norm(x)); 
    cout << "x = " << x << ", beta = " << beta << "\n";
    MTL_THROW_IF(!close(beta, 17.1464), mtl::runtime_error("wrong two_norm"));
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_infinity_norm(x)); 
    cout << "x = " << x << ", beta = " << beta << "\n";
    MTL_THROW_IF(!close(beta, 7), mtl::runtime_error("wrong one_norm"));
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_sum(x)); 
    cout << "x = " << x << ", beta = " << beta << "\n";
    MTL_THROW_IF(!close(beta, 42), mtl::runtime_error("wrong sum"));
    
    (lazy(x)= 7.0) || (lazy(beta)= lazy_product(x)); 
    cout << "x = " << x << ", beta = " << beta << "\n";
    MTL_THROW_IF(!close(beta, 117649), mtl::runtime_error("wrong sum"));
    
    (lazy(x)= 2.0) || (lazy(gamma)= lazy_dot(r, x)); 
    cout << "x = " << x << ", gamma = " << gamma << "\n";
    MTL_THROW_IF(!close(gamma, -115.2), mtl::runtime_error("wrong dot"));
    
    (lazy(r)= alpha * q) || (lazy(rho)= lazy_dot(r, q)); 
    cout << "r = " << r << ", rho = " << rho << "\n";
    MTL_THROW_IF(!close(rho, 187.2), mtl::runtime_error("wrong dot"));

    (lazy(r)= alpha * q) || (lazy(v)= 8.6 * q) || (lazy(x)= 2.2 * q); 
    cout << "r = " << r << ", v = " << v << ", x = " << x << "\n";
    MTL_THROW_IF(!close(r[0], 15.6) || !close(v[0], 17.2) || !close(x[0], 4.4), 
		 mtl::runtime_error("wrong vector scaling"));

    return 0;
}
