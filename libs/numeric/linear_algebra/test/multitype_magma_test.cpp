#include <vector>
#include <iostream>


template <typename E1, typename E2>
struct add_expr
{
    add_expr(E1 const& e1, E2 const& e2) : e1(e1), e2(e2) {}

    // std::Addable<E1, E2>::result_type
    double
    operator[] (int i) const
    {
	return e1[i] + e2[i];
    }
		
    E1 const& e1;
    E2 const& e2;
};

template <typename Value>
struct my_vector
    : public std::vector<Value>
{
    typedef std::vector<Value> base;
    typedef my_vector          self;

    my_vector(std::size_t n) : base(n) {}
    my_vector(std::size_t n, Value v) : base(n, v) {}

    using base::operator=;

    template <typename E1, typename E2>
    self& operator= (add_expr<E1, E2>const& expr)
    {
	for (int i= 0; i < this->size(); i++)
	    (*this)[i]= expr[i];
    }

    template <typename Value2>
    add_expr<self, my_vector<Value2> >
    operator+(my_vector<Value2> const& v2) const
    {
	return add_expr<self, my_vector<Value2> >(*this, v2);
    }



};

#if 0
template <typename Value1, typename Value2>
add_expr<my_vector<Value1>, my_vector<Value2> >
operator+(my_vector<Value1> const& v1, my_vector<Value2> const& v2)
{
    return add_expr<my_vector<Value1>, my_vector<Value2> >(v1, v2);
}
#endif

template <typename Value1, typename E1, typename E2>
add_expr<my_vector<Value1>, add_expr<E1, E2> >
operator+(my_vector<Value1> const& v1, add_expr<E1, E2> const& exp2)
{
    return add_expr<my_vector<Value1>, add_expr<E1, E2> >(v1, exp2);
}


template <typename Value2, typename E1, typename E2>
add_expr<add_expr<E1, E2>, my_vector<Value2> >
operator+(add_expr<E1, E2> const& expr, my_vector<Value2> const& v)
{
    return add_expr<add_expr<E1, E2>, my_vector<Value2> >(expr, v);
}

int main() 
{

    my_vector<double> v1(3, 7.0), v2(3, 8.0), v3(3, 8.0), v4(3, 9.0), v5(3, 10.0);

    v1= v2 + v3 + v4 + v4 + v5;
    std::cout << v1[0] << '\n';

    return 0;
}
