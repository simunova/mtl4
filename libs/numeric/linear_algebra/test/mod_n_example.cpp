#include <iostream>

#include <libs/numeric/linear_algebra/test/mod_n.hpp>
#include <libs/numeric/linear_algebra/test/algebraic_functions.hpp>


int main(int, char* [])
{
    using namespace mtl;
    using namespace std;

    typedef mod_n_t<unsigned, 5>    mod_5;
    math::mult<mod_5>               mult_mod_5;

    cout << "equal_results(mod_5(2), mod_5(3), mod_5(4), mod_5(4), mult_mod_5) "
	 << equal_results(mod_5(2u), mod_5(3u), mod_5(4u), mod_5(4u), mult_mod_5) << endl;

    cout << "algebraic_division(mod_5(4), mod_5(2), mult_mod_5) "
	 << algebraic_division(mod_5(4u), mod_5(2u), mult_mod_5) << endl; 

    typedef mod_n_t<unsigned, 28>    mod_28;
    math::mult<mod_28>               mult_mod_28;
    math::add<mod_28>                add_mod_28;
    
    cout << "1/3 " << math::inverse(mult_mod_28, mod_28(3)) 
	 << " check " << math::inverse(mult_mod_28, mod_28(3)) * mod_28(3) << endl;
    cout << "1/5 " << math::inverse(mult_mod_28, mod_28(5)) 
	 << " check " << math::inverse(mult_mod_28, mod_28(5)) * mod_28(5) << endl;
    cout << "1/9 " << math::inverse(mult_mod_28, mod_28(9)) 
	 << " check " << math::inverse(mult_mod_28, mod_28(9)) * mod_28(9) << endl;

    cout << "gcd(24, 28): " << gcd(24, 28) << endl;
    cout << "gcd(25, 28): " << gcd(25, 28) << endl;
    

    typedef mod_n_t<unsigned, 127>  mod_127;
    math::mult<mod_127>             mult_mod_127;
    math::add<mod_127>              add_mod_127;

    mod_127   v78(78), v113(113), v90(90), v80(80);
   
    cout << "equal_results(v78, v113,  v90, v80, mult_mod_127) "
         << equal_results(v78, v113,  v90, v80, mult_mod_127) << endl;
    cout << "equal_results(v78, v113,  v90, mod_127(81), mult_mod_127) "
         << equal_results(v78, v113,  v90, mod_127(81), mult_mod_127) << endl;

    cout << "identity_pair(v78, mod_127(-78), mult_mod_127) "
         << identity_pair(v78, mod_127(-78), mult_mod_127) << endl;
    cout << "identity_pair(v78, mod_127(57), mult_mod_127) "
         << identity_pair(v78, mod_127(57), mult_mod_127) << endl;

    cout << "algebraic_division(mod_127(8), mod_127(2), mult_mod_127) = log_2 (8) "
         << algebraic_division(mod_127(8), mod_127(2), mult_mod_127) << endl;
    cout << "algebraic_division(mod_127(35), v78, mult_mod_127) = log_78 (35) "
         << algebraic_division(mod_127(35), v78, mult_mod_127) << endl;
    
    cout << "multiply_and_square(v78, 8, mult_mod_127) = 78^8 "
	 << multiply_and_square(v78, 8, mult_mod_127) << endl;

    return 0;
}














#if 0
concept_map AbelianGroup<math::add<mod_n_t<unsigned int, 127u> >, mod_n_t<unsigned int, 127u> > {}
concept_map CommutativeSemiGroup<math::mult<mod_n_t<unsigned int, 127u> >, mod_n_t<unsigned int, 127u> > {}

concept_map GenericField
       < add< mod_n_t<unsigned, 127> >, 
	 mult< mod_n_t<unsigned, 127> >, 
	 mod_n_t<unsigned, 127> 
       > 
{
    // Why do we need the typedefs???
    
    typedef mod_n_t<unsigned, 127> inverse_result_type;
    typedef mod_n_t<unsigned, 127> identity_result_type;
    typedef bool          is_invertible_result_type;
}
#endif 


