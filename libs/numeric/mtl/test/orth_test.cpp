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
#include <cmath>
#include <complex>
#include <vector>

#include <boost/numeric/mtl/mtl.hpp>



const unsigned sz= 5;

inline float f(float x) { return x; }
inline double f(double x) { return x; }

inline std::complex<double> f(std::complex<double> x) 
{ 
    return std::complex<double>(real(x), real(x)+1.0); 
}


template <typename Vector, typename T>
void test(Vector& v, const T&, const char* name)
{
    using std::abs; using std::cout; using mtl::size; 
    using mtl::orth; using mtl::orthogonalize_factors;

    std::cout << "\n" << name << "\n";
    for (unsigned i= 0, c= 1; i < size(v); ++i)
	for (unsigned j= 0; j < size(v[i]); ++j, c++)
	    v[i][j]= f(T((i + j) % sz));


    cout << "w initially\n";
    Vector w(v);
    for (unsigned i= 0; i < size(w); ++i)
	std::cout << w[i] << "\n";
    std::cout << "\n";

    orth(w);

    for (unsigned i= 0; i < size(w); ++i)
	std::cout << w[i] << "\n";
    std::cout << "\n";

    for (unsigned i= 0, c= 1; i < size(w); ++i) {
	for (unsigned j= 0; j < size(w); ++j, ++c)
	    std::cout << dot(w[i], w[j]) << " ";
	std::cout << "\n";
    }   

    MTL_THROW_IF(abs(dot(w[3], w[4])) > 0.00001, mtl::runtime_error("Vectors 3 and 4 are not orthogonal!"));
    MTL_THROW_IF(abs(dot(w[4], w[4]) - T(1)) > 0.00001, mtl::runtime_error("Vector 4 is not normal!"));

    cout << "\nv initially\n";
    for (unsigned i= 0; i < size(v); ++i)
	std::cout << v[i] << "\n";
    std::cout << "\n";

    std::cout << "The according factors are: \n" << orthogonalize_factors(v) << '\n';

    for (unsigned i= 0; i < size(v); ++i)
	std::cout << v[i] << "\n";
    std::cout << "\n";

    for (unsigned i= 0, c= 1; i < size(v); ++i) {
	for (unsigned j= 0; j < size(v); ++j, ++c)
	    std::cout << dot(v[i], v[j]) << " ";
	std::cout << "\n";
    }   

    MTL_THROW_IF(abs(dot(v[3], v[4])) > 0.00001, mtl::runtime_error("Vectors 3 and 4 are not orthogonal!"));
    MTL_THROW_IF(abs(dot(v[4], v[4])) < 0.00001, mtl::runtime_error("Vector 4 should be non-zero!"));


}



int main(int, char**)
{
    using namespace mtl;
    dense_vector<float>                                                 cf(sz, 1.0);
    dense_vector<double>                                                cd(sz, 1.0);
    dense_vector<std::complex<double> >                                 cc(sz, 1.0);
    dense_vector<float, mtl::vec::parameters<row_major> >                 rf(sz, 1.0);

    std::vector<dense_vector<float> >                                   scf(sz, cf);
    std::vector<dense_vector<double> >                                  scd(sz, cd);
    std::vector<dense_vector<std::complex<double> > >                   scc(sz, cc);
    std::vector<dense_vector<float, mtl::vec::parameters<row_major> > >   srf(sz, rf);

    dense_vector<dense_vector<float> >                                  ccf(sz, cf);
    dense_vector<dense_vector<float>, mtl::vec::parameters<row_major> >   rcf(sz, cf);

    test(scf, cf[0], "std::vector<dense_vector<float> >");
    test(scd, cd[0], "std::vector<dense_vector<double> >");
    test(scc, cc[0], "std::vector<dense_vector<std::complex<double> > >");
    test(srf, rf[0], "std::vector<dense_vector<float, parameters<row_major> > >");

    test(ccf, cf[0], "dense_vector<dense_vector<float> >");
    test(rcf, cf[0], "dense_vector<dense_vector<float>, parameters<row_major> >");

    return 0;
}
