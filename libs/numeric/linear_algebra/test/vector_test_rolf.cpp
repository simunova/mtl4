#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>

// #include "vector_concepts.hpp"
#include <boost/numeric/linear_algebra/vector_concepts.hpp>

typedef double Type;

namespace ublas = boost::numeric::ublas;

typedef ublas::vector<Type> Vector;

namespace math {

concept_map AdditiveAbelianGroup<Vector> {
  typedef boost::numeric::ublas::vector_binary<boost::numeric::ublas::vector<double, boost::numeric::ublas::unbounded_array<double, std::allocator<double> > >, boost::numeric::ublas::vector<double, boost::numeric::ublas::unbounded_array<double, std::allocator<double> > >, boost::numeric::ublas::scalar_minus<double, double> > result_type; 
}

concept_map math::VectorSpace<Vector,Type> 
{
    typedef AdditiveAbelianGroup<Vector>::result_type          result_type;
    typedef AdditiveAbelianGroup<Vector>::assign_result_type   assign_result_type;  
}
}

template< typename Vec, typename Scalar>
requires math::VectorSpace <Vec,Scalar>
void cg(Vec& u, Scalar s) 
{

}

int main() {
  Vector u(2);
  u(0)=1.; u(1)=2.;;
  Type s=1.; 
  cg(u,s);
}
