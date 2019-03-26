/* compile with g++ ... -l arprec */

#include <iostream>
#include <fstream>
#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include <complex>
#include <boost/numeric/mtl/mtl.hpp>
#include <cassert>
#include <boost/numeric/linear_algebra/identity.hpp>

using namespace std;
typedef mtl::dense2D<mp_complex> mat_c;

ostream& operator<< (ostream& os, const mp_complex& c) {
	return os << "(" << c.real << ", " << c.imag << ")";
}

int main() {

	mp::mp_init(102);
	cout.precision(30);

	mp_real x1,x2,x3,x4,x5,x6;
	x1 = mp_real("2e18");
	x2 = mp_real("-1.34");
	x3 = mp_real("1.6");
	x4 = mp_real(3);
	x5 = mp_real("1.83e7");
	x6 = mp_real("-4.234");
	cout << "x1 = " << x1 << "\n";
	cout << "x2 = " << x2 << "\n";
	cout << "x3 = " << x3 << "\n";
	cout << "x4 = " << x4 << "\n";
	cout << "x5 = " << x5 << "\n";
	cout << "x6 = " << x6 << "\n";
	cout << "---------------\n";

	mp_complex z1,z2,z3,z4,z5,z6;
	z1 = mp_complex("4e30","3e30");
	z2 = mp_complex("3","-3.5");
	z3 = z1+z2;
	z4 = sqrt(mp_complex("-4","0"));
	z5 = z2+z4;
	z6 = z2-z3+z1+z4;
	cout << "z1 = " << z1 << "\n";
	cout << "z2 = " << z2 << "\n";
	cout << "z3 = " << z3 << "\n";
	cout << "z4 = " << z4 << "\n";
	cout << "z5 = " << z5 << "\n";
	cout << "z6 = " << z6 << "\n";
	cout << "---------------\n";
	{
	mtl::dense2D<mp_real> A(2,2),B(2,2),C(2,2),D(2,2);
	mtl::dense_vector<mp_real> b(2),x(2),r(2);
	A = x1,x2,
	    x3,x4;
	b = x5,x6;
	cout << "A = \n" << A;
	cout << "b = " << b << endl;

	// A = A*A;	// Laufzeitfehler
	A*=A;
	cout << "A = \n" << A;
	cout << "A*A = \n" << A*A;
	B =A*A;
	cout << "B = \n" << B;
	C = trans(A);
	cout << "C = \n" << C;
	D = A;

	mtl::dense_vector<int> p(2);
	lu(A,p);
	x = lu_apply(A,p,b);
	r = D*x-b;
	cout << "Ax-b = " << r << endl;
	cout << "---------------\n";
	}

	{
	mtl::dense2D<mp_complex> A(2,2),B(2,2),C(2,2),D(2,2);
	mtl::dense_vector<mp_complex> b(2),x(2),r(2);
 	A = z1,z2,	
 	    z3,z4;


	b = z5,z6;
	cout << "A = \n" << A;
	cout << "b = " << b << endl;

	// A = A*A;	// Compilezeitfehler
	// A*=A;	// Compilezeitfehler
	B =A*A;	// Compilezeitfehler
	C = trans(A);	// Compilezeitfehler

	D = A;

	mtl::dense_vector<int> p(2);
	lu(A,p);	// Compilezeitfehler
	x = lu_apply(A,p,b);
	r = D*x-b;
 	cout << "Ax-b = " << r << endl;
	cout << "---------------\n";

	cout << "abs(A[0][0]) is " << abs(A[0][0]) << '\n';
	}

	mp::mp_finalize();

}
