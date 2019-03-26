// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University. 
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG, www.simunova.com. 
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also tools/license/license.mtl.txt in the distribution.

#include <iostream>
#include <complex>

#include <boost/type_traits.hpp>

#define MTL_VERBOSE_TEST                                   //  To print out everything

// #include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/operations.hpp>
#include <boost/numeric/mtl/io/test_ostream.hpp>

// using namespace std;

using namespace mtl::io;

template <typename Vector>
struct is_float_vec
  : boost::is_floating_point<typename mtl::Collection<Vector>::value_type>
{};


template <typename Vector>
struct is_int_vec
  : boost::is_integral<typename mtl::Collection<Vector>::value_type>
{};

static int errors = 0, tests = 0;

void expect(bool cond, std::string msg)
{
    if (!cond) {
        std::cerr << "ERROR: " << msg << std::endl;
        ++errors;
    }
    ++tests;
}

template <typename Vector>
void float_only_test(Vector&, const char*, boost::integral_constant<bool, false>)
{
}

template <typename Vector>
void float_only_test(Vector& u, const char* name, boost::integral_constant<bool, true>)
{    
    Vector v(5);

    for (int i = 0; i < 5; ++i)
        u[i] = i / 4.0;
    tout << "Non-complex-" << name << ": u = " << u << "\n";    
    
    // Inverse trigonometric functions    
    v = acos(u);
    tout << "acos(u) = " << v << "\n";
    expect(std::abs(v[0] - M_PI/2.0) < 0.0001, "acos(0) should be pi/2.");
        
    v = asin(u);
    tout << "asin(u) = " << v << "\n";
    expect(std::abs(v[4] - M_PI/2.0) < 0.0001, "asin(1) should be pi/2.");
    
    v = atan(u);
    tout << "atan(u) = " << v << "\n";
    expect(std::abs(v[4] - M_PI/4.0) < 0.0001, "atan(1) should be pi/4.");
        
# ifdef MTL_WITH_MATH_ELEVEN    
    for (int i = 0; i < 5; ++i)
        u[i] = 1. + i / 2.0;
    tout << "u = " << u << "\n";
    v = acosh(u);
    tout << "acosh(u) = " << v << "\n";
    expect(std::abs(v[2] - 1.31696)  < 0.0001, "acosh(2) should be about 1.31696.");
        
    v = asinh(u);
    tout << "asinh(u) = " << v << "\n";
    expect(std::abs(v[2] - 1.44363)  < 0.0001, "asinh(2) should be about 1.44363.");
        
    for (int i = 0; i < 5; ++i)
        u[i] = -0.8 + i * 0.4;
    tout << "u = " << u << "\n";
    v = atanh(u);
    tout << "atanh(u) = " << v << "\n";
    expect(std::abs(v[4] - 1.098612289)  < 0.0001, "atanh(0.8) should be about 1.098612289.");    
# endif
    
    // Rounding operations    
    for (int i = 0; i < 5; ++i)
        u[i] = -1. + 0.4 * i;
    tout << "u = " << u << "\n";
        
    v = ceil(u);
    tout << "ceil(u) = " << v << "\n";
    expect(v[1] == 0.0, "ceil(-0.6) should be 0.");
    expect(v[4] == 1.0, "ceil(0.6) should be 1.");
    
    v = floor(u);
    tout << "floor(u) = " << v << "\n";
    expect(v[1] == -1.0, "floor(-0.6) should be -1.");
    expect(v[4] == 0.0, "floor(0.6) should be 0.");
    
# ifdef MTL_WITH_MATH_ELEVEN    
    v = round(u);
    tout << "round(u) = " << v << "\n";
    expect(v[1] == -1.0, "round(-0.6) should be -1.");
    expect(v[4] == 1.0, "round(0.6) should be 1.");
    
    v = trunc(u);
    tout << "trunc(u) = " << v << "\n";
    expect(v[1] == 0.0, "trunc(-0.6) should be 0.");
    expect(v[4] == 0.0, "trunc(0.6) should be 0.");
    
# endif
    
}

template <typename Vector>
void non_int_test(Vector&, const char*, boost::integral_constant<bool, true>)
{
}

template <typename Vector>
void non_int_test(Vector& u, const char* name, boost::integral_constant<bool, false>)
{
    Vector v(5);
    // Trigonometric functions
    for (int i = 0; i < 5; ++i)
        u[i] = i * M_PI / 8.0;
    tout << name << ": u = " << u << "\n";
    
    v = acos(u);
    tout << "cos(u) = " << v << "\n";
    expect(std::abs(v[2] - 0.667457) < 0.0001, "cos(pi/4) should be 0.667457.");
        
    v = sin(u);
    tout << "sin(u) = " << v << "\n";
    expect(std::abs(v[4] - 1.0) < 0.0001, "asin(pi/2) should be 1.");
        
    v = tan(u);
    tout << "tan(u) = " << v << "\n";
    expect(std::abs(v[2] - 1) < 0.0001, "tan(pi/4) should be 1.");
        
    v = cosh(u);
    tout << "acosh(u) = " << v << "\n";
    expect(std::abs(v[4] - 2.50918)  < 0.0001, "acosh(pi/2) should be about 2.50918.");
    
    v = sinh(u);
    tout << "sinh(u) = " << v << "\n";
    expect(std::abs(v[4] - 2.3013)  < 0.001, "sinh(pi/2) should be about 2.3013.");
    
    v = tanh(u);
    tout << "tanh(u) = " << v << "\n";
    expect(std::abs(v[4] - 0.917152)  < 0.0001, "tanh(pi/2) should be about 0.917152.");        
    
    iota(u, 1);
    tout << name << ": u = " << u << "\n";
    
    v = log(u);
    tout << "log(u) = " << v << "\n";
    expect(std::abs(v[4] - 1.60944)  < 0.0001, "log(5) should be about 1.60944.");        
    
# ifdef MTL_WITH_MATH_ELEVEN    
    v = log10(u);
    tout << "log10(u) = " << v << "\n";
    expect(std::abs(v[4] - 0.69897)  < 0.0001, "log10(5) should be about 0.69897.");        
# endif
    
}  

template <typename Vector>
void int_float_test(Vector& u, const char* name, boost::integral_constant<bool, true>)
{
    iota(u, -2);
    Vector v(5);
    tout << name << ": u = " << u << "\n";
    
    // Test that integer values should invariant regarding rounding operations
    v = ceil(u);
    tout << "ceil(u) = " << v << "\n";
    expect(v[1] == -1.0, "ceil(-2) should be -1.");
    expect(v[4] == 2.0, "ceil(2) should be 2.");
    
    v = floor(u);
    tout << "floor(u) = " << v << "\n";
    expect(v[1] == -1.0, "floor(-2) should be -1.");
    expect(v[4] == 2.0, "floor(2) should be 2.");
    
# ifdef MTL_WITH_MATH_ELEVEN    
    v = round(u);
    tout << "round(u) = " << v << "\n";
    expect(v[1] == -1.0, "round(-2) should be -1.");
    expect(v[4] == 2.0, "round(2) should be 2.");     
    
    v = trunc(u);
    tout << "trunc(u) = " << v << "\n";
    expect(v[1] == -1.0, "trunc(-2) should be -1.");
    expect(v[4] == 2.0, "trunc(2) should be 2.");     
    
    iota(u, 1);
    
    v = log2(u);
    tout << "log2(u) = " << v << "\n";    
    typename Vector::value_type log2_ref = is_int_vec<Vector>::value  ? 2.0 : 2.32193;
    expect(std::abs(v[4] - log2_ref)  < 0.0001, "log2(5) wrong."); // not elegant
    
# endif

}

template <typename Vector>
void int_float_test(Vector&, const char*, boost::integral_constant<bool, false>)
{
}


template <typename Vector>
void all_test(Vector& u, const char* name)
{
    iota(u, -2);
    Vector v(5);

    tout << name << ": u = " << u << "\n";
    v = abs(u);
    tout << "abs(u) = " << v << "\n";
    expect(v[0] == 2.0, "abs(-2) should be 2");
    expect(v[4] == 2.0, "abs(2) should be 2");
}

    
template <typename Vector>
void test(const char* name)
{
    Vector u(5);
    all_test(u, name);
    int_float_test(u, name, boost::integral_constant<bool, is_int_vec<Vector>::value || is_float_vec<Vector>::value>());
    non_int_test(u, name, is_int_vec<Vector>());
    float_only_test(u, name, is_float_vec<Vector>());
    
    tout << "\n\n";
}
 
 
 
 
int main(int, char**)
{
    using mtl::vec::parameters;
    using namespace mtl;

    tout << "Testing vector operations\n\n";

    test<dense_vector<int> >("test int");
    test<dense_vector<float> >("test float");
    test<dense_vector<double> >("test double");
# ifdef MTL_WITH_MATH_ELEVEN    
    // acos is not defined for complex in C++98
    test<dense_vector<std::complex<double> > >("test complex<double>");
# endif
    test<dense_vector<float, parameters<row_major> > >("test float in row vector");
    
    std::cout << errors << " encountered in " << tests << ".\n";
    return errors == 0 ? 0 : 1;
}
