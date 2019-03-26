#include <iostream>
#include <cmath>
#include <set>


#ifdef __GXX_CONCEPTS__
#  include <concepts>
#  include <boost/numeric/linear_algebra/new_concepts.hpp>
#else 
#  include <boost/numeric/linear_algebra/pseudo_concept.hpp>
#endif



struct mat {};                 // Matrix type
struct smat : public mat {};   // Symmetric matrix type


concept Symmetric<typename Matrix> { /* axioms */ }

concept_map Symmetric<smat> {}

std::set<const void*> table_Symmetric;

template <typename Matrix>
void inline map_Symmetric(const Matrix& A) 
{ table_Symmetric.insert(&A); }

template <Symmetric Matrix>
void inline map_Symmetric(const Matrix&) {}

template <typename Matrix>
void inline unmap_Symmetric(const Matrix& A) 
{ table_Symmetric.erase(&A); }

template <Symmetric Matrix>
void inline unmap_Symmetric(const Matrix&) {}

template <typename Matrix>
bool inline is_Symmetric(const Matrix& A) 
{ return table_Symmetric.find(&A) != table_Symmetric.end(); }

template <Symmetric Matrix>
bool inline is_Symmetric(const Matrix&) 
{ return true; }




concept PositiveDefinite<typename Matrix>{ /* axioms */ }

std::set<const void*> table_PositiveDefinite;

template <typename Matrix>
void inline map_PositiveDefinite(const Matrix& A) 
{ table_PositiveDefinite.insert(&A); }


template <typename Matrix>
void inline unmap_PositiveDefinite(const Matrix& A) 
{ table_PositiveDefinite.erase(&A); }

template <typename Matrix>
bool inline is_PositiveDefinite(const Matrix& A) 
{ return table_PositiveDefinite.find(&A) != table_PositiveDefinite.end(); }

template <PositiveDefinite Matrix>
bool inline is_PositiveDefinite(const Matrix&) 
{ return true; }






template <typename Matrix>
void spd_solver(const Matrix& A)
{
    std::cout << "spd_solver (Symmetric positiv-definit)\n";
}

template <typename Matrix>
void symmetric_solver(const Matrix& A)
{
    std::cout << "symmetric_solver\n";
}

template <typename Matrix>
void default_solver(const Matrix& A)
{
    std::cout << "Default_solver\n";
}

template <typename Matrix>
void solver(const Matrix& A)
{
    if (is_Symmetric(A)) {
	if (is_PositiveDefinite(A))
	    spd_solver(A);
	else
	    symmetric_solver(A);
	return;
    }

    default_solver(A);
}

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
