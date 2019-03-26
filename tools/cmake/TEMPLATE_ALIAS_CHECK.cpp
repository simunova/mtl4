// #include <iostream>


template <typename T, typename U>
struct dings {};

template <typename T>
using bums= dings<T, T>;

int main()
{
    typedef bums<int> bi;

    return 0;
}
