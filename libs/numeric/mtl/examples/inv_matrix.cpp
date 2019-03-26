#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main(int, char**)
{
    using namespace mtl; using namespace std;
    typedef  mtl::dense2D<double>            Matrix;

    double array[][4]= {{2, 3,   4,   5}, 
                        {4, 10, 13,  16},
                        {6, 25, 38,  46},
		        {8, 32, 77, 100}};
    Matrix 		A(array), I(4, 4);
    I= 1.0;

    Matrix LU(A);
    dense_vector<std::size_t> v(4);
    lu(LU, v);
    mat::traits::permutation<>::type P(permutation(v));
    
    cout << "A is:\n" << A << "\nPermuted A is \n" << Matrix(P * A);

    Matrix L(I + strict_lower(LU)), U(upper(LU)), A2(L * U);
    cout << "L [permuted] is:\n" << L << "U [permuted] is:\n" << U 
	 << "L * U [permuted] is:\n" << A2
	 << "L * U is:\n" << Matrix(trans(P) * A2);
 
    Matrix UI(inv_upper(U));
    cout << "inv(U) [permuted] is:\n" << UI << "UI * U is:\n" << UI * U;
 
    Matrix LI(inv_lower(L));
    cout << "inv(L) [permuted] is:\n" << LI << "LI * L is:\n" << LI * L;
 
    Matrix AI(UI * LI * P);
    cout << "inv(A) [inv(U) * inv(L) * P] is \n" << AI << "Test: A * AI is\n" << AI * A;
 
    mat::traits::inv<Matrix>::type A_inv(inv(A));
    cout << "inv(A) is \n" << A_inv << "Test: A * AI is\n" << A_inv * A;

    return 0;
}
