#include <iostream>

#include <boost/numeric/linear_algebra/operators.hpp>
#include <boost/numeric/linear_algebra/concepts.hpp>

#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>

int main(int, char* []) 
{
    using namespace mtl;

    math::add<int>     int_add;

    std::cout << "equal_results(2,5,  3,4, int_add) " 
         << equal_results(2,5,  3,4, int_add)  << std::endl;
    std::cout << "equal_results(2,4,  3,4, int_add) " 
         << equal_results(2,4,  3,4, int_add)  << std::endl;

    std::cout << "identity_pair(2,4, int_add) " 
         << identity_pair(2,4, int_add)  << std::endl;
    std::cout << "identity_pair(0,0, int_add) " 
         << identity_pair(0,0, int_add)  << std::endl << std::endl;

    std::cout << "ordered_algebraic_division(3,2, int_add) " 
         << ordered_algebraic_division(3,2, int_add)  << std::endl;
    std::cout << "ordered_algebraic_division(4,2, int_add) " 
         << ordered_algebraic_division(4,2, int_add)  << std::endl;
    std::cout << "ordered_algebraic_division(-6,2, int_add) " 
         << ordered_algebraic_division(-6,2, int_add)  << std::endl;

    return 0;
}
 
