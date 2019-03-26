#include <iostream>                                                                                             
#include <typeinfo>                                                                                             
#include <boost/numeric/mtl/mtl.hpp>                                                                            
#include <boost/type_traits/is_const.hpp>                                                                            
                                                                                                                
using namespace mtl; using mtl::iall;                                                                           
                                                                                                                
template <typename T>                                                                                                 
void beast2(T& B)
{
    B[1][1] = 666;
}

//                                                                                                              
// we feed the beast with a const reference                                                                     
//                                                                                                              
void                                                                                                            
beast(const dense2D<double> &A)                                                                                 
{                                  
    // sub_matrix(A, 0, 5, 0, 5)= "Hallo";
    auto B = sub_matrix(A, 0, 5, 0, 5);                             
    // dense2D<double>   B = A[irange(5)][irange(5)];
    B[1][1] = 666;            
    // sub_matrix(A, 0, 5, 0, 5)[1][1] = 666; // const-aware error message

    beast2(B);
    // beast2(sub_matrix(A, 0, 5, 0, 5));  // const-aware error message
}         
     
template <typename T>
void print(T&)
{
    std::cout << boost::is_const<T>::value << '\n';
}

template <typename T>
struct bums {};

template <typename T>
bums<T> make_bums(T&)
{
    return bums<T>();
}
    
class dings
{
  public:
    dings()
    {
	std::cout << typeid(*this).name() << '\n';
	std::cout << boost::is_const<decltype(*this)>::value << '\n';
	print(*this);
	std::cout << typeid(make_bums(*this)).name() << '\n';
    }

    void quatsch() const { print(*this); }
    void quatsch()       { print(*this); }
};

                                                                                                            
int main(int, char**)                                                                                           
{ 
    dings       a;
    const dings b;

    std::cout << boost::is_const<decltype(b)>::value << '\n';

    print(a);
    print(b);

    std::cout << typeid(a).name() << '\n';
    std::cout << typeid(b).name() << '\n';

    std::cout << typeid(make_bums(a)).name() << '\n';
    std::cout << typeid(make_bums(b)).name() << '\n';

    a.quatsch();
    b.quatsch();

    return 0;                                                                                                              
    using namespace mtl; using mtl::iall;                                                                       
                                                                                                                
    const int         m = 5, n = 5;                                                                             
    dense2D<double>   A(n, n);                                                                                  
                                                                                                                
    for (int i=0; i<m; ++i) {                                                                                   
        for (int j=0; j<n; ++j) {                                                                               
            A[i][j] = i+j;                                                                                      
        }                                                                                                       
    }                                                                                                           
    const dense2D<double>   B(A);                                                 
    dense2D<double>   C(B);
                                               
    std::cout << "before: A is\n" << A << "\n";                                                                 
//                                                                                                              
//  calling the beast                                                                                           
//                                                                                                              
    beast(A);          
    std::cout << "after: A is\n" << A << "\n";                                                                  
                       
                                                                                         
    return 0;                                                                                                   
}
