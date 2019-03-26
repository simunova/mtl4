// Filename: qr_givens_example.cpp (part of MTL4)

#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

using namespace std;
typedef mtl::mat::dense2D<double> dMatrix;

int main() {
    dMatrix M1(3,3), M2(3,3);
    

    M1 = 2,0,0,
	1,1,0,
	0,1,3; //EWs: 1,2,3     
	
    mtl::mat::qr_givens_solver<dMatrix> QR1(M1);
    QR1.setTolerance(1.0e-5);
    QR1.calc();
    cout << "M1(providing tolerance):\n Q: \n"  << QR1.getQ() << "\n R: \n" << QR1.getR() << "\n";      
    M2 = -261, 209,  -49,
        -530, 422,  -98,
        -800, 631, -144; //EWs: 3,4,10  
        
    mtl::mat::qr_givens_solver<dMatrix> QR2(M2);
    QR2.calc();
    cout << "M2(with defaults):\n Q: \n"  << QR2.getQ() << "\n R: \n" << QR2.getR() << "\n";
    
    dMatrix Q2, R2;
    boost::tie(Q2, R2)= qr_givens(M2);
    cout << "M2(with defaults):\n Q: \n"  << Q2 << "\n R: \n" << R2 << "\n";

    return 0;
}
