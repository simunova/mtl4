#include <iostream>
#include <cmath>


concept C1<typename T> {}
concept C2<typename T> {}
concept C3<typename T> {}

template <typename T>
  requires C3<T>
concept_map C2<T> {}

template <typename T>
  requires C1<T>
concept_map C2<T> {}

template <typename T>
  requires C2<T>
concept_map C1<T> {}

concept_map C1<int> {}
concept_map C2<short> {}
concept_map C3<char> {}


template <C1 T> 
void h(const T& x) {}

int main(int, char* [])  
{
    int   i;
    short s;
    char  c;
    long  l;

    h(i);
    h(s);
    h(c);
    // h(l); hoped for infinite loop

    return 0;
}

