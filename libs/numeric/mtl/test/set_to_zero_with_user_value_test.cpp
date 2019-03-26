#include <boost/numeric/mtl/mtl.hpp>

template < class T >
struct MyContainer
{
    MyContainer( size_t dim = 0 ) : dim_(dim) {}

    size_t dim_;
};

template < class T >
MyContainer<T> operator * ( MyContainer<T> const& a, MyContainer<T> const& )
{
    return MyContainer<T>(a.dim_);
};

template < class T >
MyContainer<T>& operator += ( MyContainer<T>& a, MyContainer<T> const&)
{
    return a;
};

template < class T >
void set_to_zero( MyContainer<T>&) {
    // ... //
}

namespace mtl {
namespace ashape {

template < class T >
struct ashape< MyContainer<T> > : nonscal {
    typedef nonscal type;
};

} // namespace ashape
} // namespace mtl


int main( void ) 
{
    typedef MyContainer<double> Container;
    typedef mtl::dense2D<Container> Matrix;

    Matrix m1(2,2);

    m1(0,0) = Container(2);
    m1(0,1) = Container(2);
    m1(1,0) = Container(2);
    m1(1,1) = Container(2);

    std::cout << "m1 dim = " << m1(0,0).dim_ << std::endl;

    Matrix m2(2,2);

    m2(0,0) = Container(2);
    m2(0,1) = Container(2);
    m2(1,0) = Container(2);
    m2(1,1) = Container(2);

    Matrix m3(2,2);

    m3(0,0) = Container(2);
    m3(0,1) = Container(2);
    m3(1,0) = Container(2);
    m3(1,1) = Container(2);

    mult(m1,m2,m3);  // calls a tiling functor, is this default?

    // mtl::gen_dmat_dmat_mult_ft<Matrix,Matrix,Matrix>()(m1,m2,m3);

    // m3= m1 * m2; // meta-program kommt nicht mit non_scal klar

    //mtl::gen_dmat_dmat_mult_ft<Matrix,Matrix,Matrix,mtl::assign::plus_sum>()(m1,m2,m3);

    std::cout << "m3 dim = " << m3(0,0).dim_ << std::endl;

    MTL_THROW_IF(m3(0,0).dim_ != 2, mtl::runtime_error("Dimension of matrix element overwritten"));

    return 0;
}

