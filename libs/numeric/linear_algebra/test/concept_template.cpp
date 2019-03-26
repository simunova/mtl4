#include <iostream>
#include <cmath>
#include <set>


#ifdef __GXX_CONCEPTS__
#  include <concepts>
#  include <boost/numeric/linear_algebra/new_concepts.hpp>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif



template <concept C>
class dynamic_concept
{
    static std::set<void*> table;    
  public:
    template <typename T>
    bool is(const T& x) { return table.find(&x) != table.end(); }

    template <C T>
    bool is(const T& x) { return true; }

    template <typename T>
    void map(const T& x) { table.insert(&x); }
    
    template <C T> void map(const T& x) {}
    
    template <typename T>
    void un_map(const T& x) { table.erase(&x); }
	
    template <C T> void unmap(const T& x) {}
};




int main(int, char* [])  
{




    return 0;
}
