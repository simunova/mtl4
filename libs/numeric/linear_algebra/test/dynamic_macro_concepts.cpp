#define DYNAMIC_CONCEPT(SCONCEPT)		\
    std::set<const void*> table_ ## SCONCEPT;   \
                                                \
    template <typename T>                       \
    void inline map_ ## SCONCEPT(const T& x)    \
    { table_ ## SCONCEPT.insert(&x); }          \
                                                \
    template <SCONCEPT T>                       \
    void inline map_ ## SCONCEPT(const T& x) {}	\
                                                \
    template <typename T>                       \
    void inline unmap_ ## SCONCEPT(const T& x)  \
    { table_ ## SCONCEPT.erase(&x); }           \
                                                \
    template <SCONCEPT T>                       \
    void inline unmap_ ## SCONCEPT(const T& x) {}	\
                                                \
    template <typename T>                       \
    bool inline is_ ## SCONCEPT(const T& x)     \
    { return table_ ## SCONCEPT.find(&x) != table_ ## SCONCEPT.end(); } \
                                                \
    template <SCONCEPT T>                       \
    bool inline is_ ## SCONCEPT(const T&)       \
    { return true; }                            


class mat {};                 // Matrix type
class smat : public mat {};   // Symmetric matrix type


concept Symmetric<typename Matrix> { /* axioms */ }
concept PositiveDefinite<typename Matrix>{ /* axioms */ }

DYNAMIC_CONCEPT(Symmetric)
DYNAMIC_CONCEPT(PositiveDefinite)

concept_map Symmetric<smat> {}


template <typename Matrix>
void spd_solver(const Matrix& A)
{  std::cout << "SPD_solver (Symmetric and positiv-definite)\n"; }

template <typename Matrix>
void symmetric_solver(const Matrix& A)
{ std::cout << "Symmetric_solver\n"; }

template <typename Matrix>
void solver(const Matrix& A)
{
    if (is_PositiveDefinite(A) && is_Symmetric(A)) { spd_solver(A); return; }
    if (is_Symmetric(A)) { symmetric_solver(A); return; }

    std::cout << "Default_solver\n";
}



#if 0 // What I'm emulating here:
      // This can be sorted in the concept lattice
      // and it can deal with different arities and return types

template <typename Matrix>
void solver(const Matrix& x)
{    std::cout << "Default_solver\n"; }

template <typename Matrix>
void solver(const Matrix& A)
    requires Symmetric(A);
{ std::cout << "Symmetric_solver\n"; }

template <typename Matrix>
void solver(const Matrix& A)
    requires Symmetric(A) 
          && PositiveDefinite(A);
{  std::cout << "SPD_solver (Symmetric and positiv-definit)\n"; }

#endif 


int main(int, char* [])  
{
    mat  A, B;
    smat C, D;

    map_Symmetric(B);
    map_PositiveDefinite(D);

    solver(A);
    solver(B);
    solver(C);
    solver(D);

    return 0;
}







#if 0 // The output (as expected):

Default_solver
Symmetric_solver
Symmetric_solver
SPD_solver (Symmetric positiv-definite)

#endif
