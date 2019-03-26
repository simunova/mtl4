#include <iostream>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>

#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>

// float has a concept map for Field 
// and PartiallyInvertibleMonoid w.r.t. mult
// This implies concept maps for all other scalar concepts

using math::identity; using math::add;

int main(int, char* [])
{
    using namespace std;
    using namespace mtl;

    math::add<float>    float_add;
    math::mult<float>   float_mult;

    cout << "equal_results(2.,5.,  3.,4., float_add) " 
         << equal_results(2.f,5.f,  3.f,4.f, float_add)  << endl;
    cout << "equal_results(2.,6.,  3.,4., float_add) " 
         << equal_results(2.,6.,  3.,4., float_add)  << endl;

    cout << "identity_pair(2.,4., float_mult) "  
         << identity_pair(2.f,4.f, float_mult)  << endl;
    cout << "identity_pair(0.5,2., float_mult) " 
         << identity_pair(0.5f,2.f, float_mult)  << endl << endl;

    cout << "ordered_algebraic_division(32.,2., float_mult) " 
         << ordered_algebraic_division(32.f,2.f, float_mult)  << endl;
    cout << "ordered_algebraic_division(33.f,2.f, float_mult) " 
         << ordered_algebraic_division(33.f,2.f, float_mult)  << endl;
    cout << "ordered_algebraic_division(33.f,2.f, float_add) " 
         << ordered_algebraic_division(33.f,2.f, float_add)  << endl;
    cout << "ordered_algebraic_division(31.,2., float_mult) " 
         << ordered_algebraic_division(31.f,2.f, float_mult)  << endl;
    cout << "ordered_algebraic_division(0.125,2, float_mult) " 
         << ordered_algebraic_division(0.125f,2.f, float_mult)  << endl;
    
  return 0;
}
