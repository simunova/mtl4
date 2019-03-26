#include <cstddef>

// Should be extended to be sure that support is sufficient

template <typename ValueType, ValueType FirstValue, ValueType ...Values>
struct static_vector
{
    static const std::size_t length= sizeof...(Values) + 1;
};

template <typename ValueType, ValueType FirstValue>
struct static_vector<ValueType, FirstValue>
{};

int main()
{
    typedef static_vector<int, 3, 4> sv;

    return 0;
}
