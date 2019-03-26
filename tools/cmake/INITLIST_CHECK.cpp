#include <initializer_list>

struct ctor {};

template <typename V>
struct my_vector
{
    typedef my_vector<V> self;

    template <typename T>
    my_vector(std::initializer_list<T> ls) 
    {
	for (auto l : ls);
    }

    my_vector(self&&) {}
    my_vector(const self&) {}
    my_vector(const self&, ctor) {}

    self& operator=(self&&) {}
};

int main() 
{
    my_vector<double> v= {3, 4, 5};
    v= {4, 5, 6};

    const my_vector<double> w= {3, 4, 5};

    return 0;
}
