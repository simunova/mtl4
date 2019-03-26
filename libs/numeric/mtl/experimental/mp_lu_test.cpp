/* compile with g++ mp_test.cpp -o mp_test -l arprec */
// g++ mp_lu_test.cpp -o mp_lu_test -I/home/pgottsch/Software/arprec-2.2.7/include -I$BOOST_ROOT -I$MTL -L/home/pgottsch/Software/arprec-2.2.7/src -l arprec -DMTL_HAS_ARPREC

#include <iostream>
#include <arprec/mp_real.h>
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;

int main() {

	mp::mp_init(102);

	mp_real a,b,c;

	a = "1.0e30";
	b = a+5;
	c = b-mp_real("1e30");
	c/= b;

	cout.precision(50);
	cout << "a = " << a << "\n";
	cout << "b = " << b << "\n";
	cout << "c = " << c << "\n";
	cout << "---------------\n";

	double d,e,f;
	d = 1e30;
	e = d+5;
	f = e-1e30;
	cout << "d = " << d << "\n";
	cout << "e = " << e << "\n";
	cout << "f = " << f << "\n";
	cout << "---------------\n";

	complex<double> z1;
	z1 = complex<double>(4,3);
	cout << "z1 = " << z1 << "\n";
	cout << "---------------\n";

	mtl::dense2D<mp_real> A(3,3);
	A = mp_real("2.43"), mp_real("2.33"), mp_real("7.18"),
	    mp_real("5.83"), mp_real("4.91"), mp_real("3.49"),
	    mp_real("8.12"), mp_real("0.45"), mp_real("6.33");
	cout << "A = \n" << A << "\n";
	mtl::dense_vector<int> p(3);
	
#if 0
	mtl::ashape::ashape<int>::type x= "";
	mtl::ashape::ashape<mp_real>::type y= "";
	mtl::traits::is_scalar<int>::type x= "";
	mtl::traits::is_scalar<mp_real>::type y= "";
#endif
	lu(A,p);
	cout << "A = \n" << A << "\n";
}

