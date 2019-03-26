#include<boost/numeric/ublas/matrix.hpp>
#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include<complex>

#include <boost/numeric/linear_algebra/concepts.hpp>
#include <boost/numeric/linear_algebra/vector_concepts.hpp>

typedef double Type;

namespace ublas = boost::numeric::ublas;

typedef ublas::vector<Type> Vector;

template <typename T>
struct dot 
{
  T operator() (const T& v, const T& w)
  {
    using std::conj;
    T tmp= 0;
    for (int i= 0; i < v.size(); i++)
      tmp+= conj(v[i]) * w[i];
    return tmp;
  }
};


namespace math {
  concept_map HilbertSpace< dot<double>, Vector>;
}


int main (int argc, char* argv[])  
{
  const int v_size= 10;

  Vector v(v_size), w(v_size), x, y(v_size-1);
  for (int i= 0; i < v_size; i++)
    v[i]= 1.0, w[i]= 2.0;
  for (int i= 0; i < v_size-1; i++) 
    y[i]= 3.0;

  w+= v;
  try {
    x= v + y;
  } catch (ublas::bad_argument) {
    std::cout << "caught bad argument " << std::endl;
  }


  std::cout << "v " << v << std::endl;
  std::cout << "w " << w << std::endl;
  std::cout << "x " << x << std::endl;

  return 0;
}
