#include "odd_or_even.hpp"
//#include <boost/numeric/meta_math/odd_or_even.hpp>
#include<iostream>
using namespace std;
int main(){
        const long long n=6174;
        cout<<n<<" is "<<(meta_math::odd_or_even<n>::value ? "ODD":"EVEN")<<'\n';     
        return 0;
    }
