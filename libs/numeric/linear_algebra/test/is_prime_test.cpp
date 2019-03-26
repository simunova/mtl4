#include <iostream>
#include <boost/numeric/meta_math/is_prime.hpp>

int main()
{
    const long int i= 10001;


    std::cout << i << " is " << (meta_math::is_prime<i>::value ? "" : "not ") << "prime.\n";


    return 0;
}
