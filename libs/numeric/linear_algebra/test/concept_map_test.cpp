#include <concepts>
#include <iostream>
#include <string>

#if 0 // Crashes ConceptGCC :-!
concept CExplicit<typename T>
{
    requires std::CopyConstructible<T>;
    requires std::CopyAssignable<T>;
    // How do I request printability????
    std::ostream& operator<<(std::ostream& os, T const& t);
}
template <CExplicit<T> >
void f(const T& x)
{
    std::cout << "In f with x = " << x << "\n";
}
#endif

concept CExplicit<typename T>
{
    requires std::CopyConstructible<T>;
    requires std::CopyAssignable<T>;
}


template <CExplicit T>
void f(const T& x, const char* xc)
{
    T y(x);
    std::cout << "In f with x = " << xc << "\n";
}



// Now faking concept_map
auto concept CAuto<typename T>
{
    //requires std::CopyConstructible<T>;
    //requires std::CopyAssignable<T>;
}
template <CAuto T> concept_map CExplicit<T> {}



int main() 
{
    f(3, "3");
    f(3.14, "3.14");
#if 0 // Another crash source in ConceptGCC
    f("Godmorgen. Har du sove godt?", "Godmorgen. Har du sove godt?");
#endif
    f(std::string("Jeg kan ikke klage."), "Jeg kan ikke klage.");

    return 0;
}
