#include <iostream>

#include <libs/numeric/linear_algebra/test/mod_n.hpp>
#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>
#include <libs/numeric/linear_algebra/test/power.hpp>

using mtl::mod_n_t;

template <typename Op, typename Element>
void compute_power(Element base, int exp, Op op, std::string op_name)
{
    using mtl::power;
    try {
	std::cout << op_name << ":   " << base << "^" << exp << " == " << "  " << power(base, exp, op) << '\n';
    } catch (char const* message) {
	std::cout << "\n==== Exception caught: " << message << '\n';
    }
}

template <typename Element>
void test(Element, std::string  mod_name)
{
    using std::string; using std::cout;

    math::add<Element>         add;
    string                     add_name(string("Add modulo ") + mod_name);

    for (int i= 13; i < 19; i++)
	compute_power(Element(i), 213, add, add_name);

    cout << '\n';
    for (int i= 13; i < 19; i++)
	compute_power(Element(i), -213, add, add_name);

    math::mult<Element>        mult;
    string                     mult_name(string("Mult modulo ") + mod_name);
     
    cout << '\n';
    for (int i= 13; i < 19; i++)
	compute_power(Element(i), 213, mult, mult_name);
 
    cout << '\n';
    for (int i= 13; i < 19; i++)
	compute_power(Element(i), -213, mult, mult_name);

    cout << '\n';
}



int main(int, char* []) 
{
    test(mod_n_t<unsigned, 16>(), "16");
    test(mod_n_t<unsigned, 17>(), "17");
    test(mod_n_t<unsigned, 18>(), "18");
    test(mod_n_t<unsigned, 19>(), "19");

    return 0;
}


